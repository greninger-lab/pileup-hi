use crate::{
    bamio::{BamDataSource, BamReader},
    baq::realign_record,
    cigar_resolve::resolve_cigar,
    output::{OrderedPileupOutput, OutputMethod},
    params::PileupParams,
    position_queue::GenomeInterval,
    read_buf::{BufPushResult, ReadBuffer},
    read_filter::ReadFilter,
    refseq::RefSeq,
    utils::read_ends_before_pos,
};

use anyhow::Error;
use rust_htslib::bam::Record;

pub struct PileupIterator<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    cur_interval: usize,
    tid: i32,
    next_tid: i32,
    last_tid_with_cov: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,

    rbuf: ReadBuffer,
    output: Option<T>,
    dest: OutputMethod<T>,
    pub reader: BamReader,
    refseq: Option<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    show_empty_coords: bool,
    show_empty_regions: bool,
    realign: bool,
    min_baseq: u8,
    min_mapq: u8,
    redo_baq: bool,
}

impl<T: OrderedPileupOutput> PileupIterator<T> {
    pub fn new(
        src: &BamDataSource,
        intervals: &[GenomeInterval],
        params: &PileupParams,
        output: T,
        dest: OutputMethod<T>,
    ) -> Result<Self, Error> {
        assert!(!intervals.is_empty());
        let reader = BamReader::new(src, 2)?;

        let rbuf = ReadBuffer::new(params.depth, params.disable_overlaps);

        let read_filter = ReadFilter::new(
            params.count_orphans,
            params.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        let cur_rec = Record::new();

        let show_empty_coords = params.show_empty_coords || params.show_empty_regions;
        let show_empty_regions = params.show_empty_regions;
        let min_baseq = params.min_baseq;
        let min_mapq = params.min_mapq;

        let refseq = if let Some(ref_file) = &params.refseq {
            Some(RefSeq::from_file(ref_file)?)
        } else {
            None
        };

        let pos @ next_pos @ max_pos = -1;
        let tid @ next_tid @ last_tid_with_cov = -1;

        Ok(Self {
            tid,
            next_tid,
            last_tid_with_cov,
            intervals: intervals.to_vec(),
            cur_interval: 0,
            pos,
            next_pos,
            max_pos,
            rbuf,
            output: Some(output),
            dest,
            reader,
            read_filter,
            refseq,
            cur_rec,
            min_baseq,
            min_mapq,
            show_empty_coords,
            show_empty_regions,
            realign: !params.no_baq,
            redo_baq: params.redo_baq,
        })
    }

    /// Generate a pileup from all bases passing the minimum quality filter and covering the
    /// iterator's current reference position.
    ///
    /// If allocate is true, allocate a new output type T
    #[inline(always)]
    pub fn set_pileup(&mut self) -> Result<(), Error> {
        assert!(self.rbuf.backup_buf.is_empty());
        let mut skip = false;

        // don't bother going through read buffer if it starts beyond the
        // current coordinate
        if !self.show_empty_regions && !self.show_empty_coords {
            if let Some((head_tid, head_pos)) = self.rbuf.head() {
                if head_tid == self.tid && head_pos > self.pos {
                    skip = true;
                }

                if head_tid > self.tid {
                    skip = true
                }
            }
        }

        let dest = &mut self.dest;

        let ref_sequence = &mut self.refseq.as_ref().and_then(|r| r.yield_seq());
        let rbuf = &mut self.rbuf;

        match dest {
            OutputMethod::WriteDirectly(ref mut writer) => {
                let mut output = self.output.take().unwrap();
                output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);

                let generated = if !skip {
                    generate_pileup(rbuf, ref_sequence, &mut output, self.pos, self.tid, self.min_baseq)?
                } else {
                    false
                };

                if generated || output.depth() > 0 || self.show_empty_coords {
                    self.last_tid_with_cov = self.tid;
                    output.write(writer)?;
                } else {
                    output.clear();
                }
                self.output = Some(output);
            }

            OutputMethod::QueueForOutput(output_chunk) => {
                let output = output_chunk.get_current_mut();
                output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);
                let generated = if !skip {
                    generate_pileup(rbuf, ref_sequence, output, self.pos, self.tid, self.min_baseq)?
                } else {
                    false
                };
                if generated || output.depth() > 0 || self.show_empty_coords {
                    self.last_tid_with_cov = self.tid;
                    output_chunk.advance();
                } else {
                    output_chunk.tombstone();
                }
            }
        }

        Ok(())
    }

    fn set_ref(&mut self, interval: GenomeInterval) -> Result<(), Error> {
        if interval.tid >= self.reader.header.target_count() as i64 {
            anyhow::bail!("Interval has TID exceeding header maximum!");
        }

        // end +1 because position-queues are 0-indexed and end fetch boundary is exclusive
        self.reader
            .init_to_ref(interval.tid as u32, interval.start, interval.end + 1)?;

        self.tid = interval.tid as i32;
        self.next_tid = self.tid;

        self.pos = interval.start;
        self.next_pos = interval.start;
        self.max_pos = interval.end;

        if let Some(refseq) = &mut self.refseq {
            refseq.load_seq(&self.reader.cur_ref)?;
        }

        Ok(())
    }

    // load the read buffer until we either 1) run out of data or 2) hit a read at the next
    // position/tid.
    #[inline(always)]
    pub fn intake(&mut self) -> Result<IterResult, Error> {
        // self.pos = self.pos.min(self.next_pos);
        // if self.pos > self.max_pos || {
        //     return Ok(IterResult::ReferenceEnd);
        // }

        while self.pos == self.next_pos || self.next_pos <= self.pos {
            // we need to keep reading until we have gathered all reads overlapping a position.
            //
            // TODO: move the IO reading logic outside
            if let Some(read) = self.reader.read_no_alloc(&mut self.cur_rec) {
                read?;
                let r = &mut self.cur_rec;

                if r.is_unmapped() {
                    continue;
                }

                if !self.read_filter.check_read(r) {
                    continue;
                }

                // we passed queried region
                if r.pos() > self.max_pos {
                    return Ok(IterResult::ReferenceEnd);
                }

                if self.realign {
                    if let Some(refseq) = self.refseq.as_ref().and_then(|r| r.yield_seq()) {
                        let flag = if self.redo_baq { 7 } else { 3 };
                        realign_record(r, refseq, refseq.len() as i64, flag)?;
                    }
                }

                if r.mapq() < self.min_mapq {
                    continue;
                }

                let ret = self.rbuf.attempt_push(self.tid, self.pos, r)?;

                match ret {
                    BufPushResult::Unmapped => panic!(),

                    BufPushResult::DifferentReference => {
                        self.next_tid = r.tid();
                        self.next_pos = r.pos();
                        return Ok(IterResult::ReferenceEnd);
                    }

                    BufPushResult::Pushed => {
                        self.next_pos = r.pos();
                    }
                    // self.next_pos = self.next_pos.max(r.pos())}

                    // if we've capped our buffer to a given depth, we'll iterate over all
                    // remaining reads spanning this coordinate before stopping to generate
                    // pileups. This way we won't have to deal with them at the next position.
                    BufPushResult::MaxDepthMet | BufPushResult::BeforePos => {
                        continue;
                    }
                }
            } else {
                // we ran out of reads.
                return Ok(IterResult::ReferenceEnd);
            }
        }

        // we're at the same TID as when we started, but we hit a read starting at the next
        // coordinate.
        Ok(IterResult::Generated)
    }

    pub fn auto_loop(&mut self) -> Result<(), Error> {
        assert!(!self.intervals.is_empty());
        self.set_ref(self.intervals[0].clone())?;

        loop {
            match self.intake()? {
                IterResult::Generated => {
                    // eprintln!("Generated {} -> {} / {}", self.pos, self.next_pos, self.max_pos);
                    while self.pos < self.next_pos {
                        self.set_pileup()?;
                        self.pos += 1;
                    }
                }

                IterResult::ReferenceEnd => {
                    // if we have reads for current ref still in buffer, process them until they no
                    // longer overlap with cur pos.
                    while let Some((head_tid, _)) = self.rbuf.head() {
                        if head_tid == self.tid && self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1;
                        } else {
                            break; // done with region
                        }
                    }

                    // if we are showing empty coords and this region has depth ANYWHERE, emit for
                    // all coords.
                    if self.last_tid_with_cov == self.tid && self.show_empty_coords {
                        while self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1;
                        }
                    }

                    // if we are showing empty regions, emit for entire ref regardless of depth.
                    if self.show_empty_regions {
                        while self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1;
                        }
                    }

                    self.cur_interval += 1;

                    if self.cur_interval >= self.intervals.len() {
                        break;
                    }

                    // if showing empty regions, repeat emission for entire refs we skipped while
                    // intaking reads.
                    if self.next_tid as i64 > self.intervals[self.cur_interval].tid + 1 && self.show_empty_regions {
                        while self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1;
                        }

                        while self.tid < self.next_tid {
                            self.set_ref(self.intervals[self.cur_interval].clone())?;

                            while self.pos <= self.max_pos {
                                self.set_pileup()?;
                                self.pos += 1;
                            }

                            self.cur_interval += 1;
                        }
                    }

                    if self.cur_interval >= self.intervals.len() {
                        break;
                    } else {
                        self.set_ref(self.intervals[self.cur_interval].clone())?;
                    }
                }
            }
        }

        // if we are storing output in intermediate buffer, flush it.
        match &mut self.dest {
            OutputMethod::WriteDirectly(_) => (),
            OutputMethod::QueueForOutput(out) => out.flush(),
        }

        Ok(())
    }
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
}

#[inline(always)]
pub fn generate_pileup<T: OrderedPileupOutput>(
    rbuf: &mut ReadBuffer,
    ref_sequence: &Option<&[u8]>,
    out: &mut T,
    pos: i64,
    tid: i32,
    min_baseq: u8,
) -> Result<bool, Error> {
    let mut skip_remainder_of_buf = false;
    let mut generated = false;

    for raw in rbuf.rbuf.drain(..) {
        // from a previous record, we decided to skip all remaining records in this buffer.
        if skip_remainder_of_buf {
            rbuf.backup_buf.push(raw);
            continue;
        }

        let mut r = raw.borrow_mut();

        // record starts beyond position, which means that the remainder of the buffer does
        // too. Skip the rest of the records.
        if r.rec.tid() > tid || (r.rec.pos() > pos && r.rec.tid() == tid) {
            // println!("DISCARDING: {} {} | {} {}", r.rec.pos(), r.rec.tid(), pos, tid);
            drop(r);
            rbuf.backup_buf.push(raw);
            skip_remainder_of_buf = true;
            continue;
        }

        // record is old and no longer overlaps the query coordinate. Discard.
        if read_ends_before_pos(&r, pos) || r.rec.tid() < tid {
            rbuf.depth -= 1;
            continue;
        }

        generated = true;

        // advance to the current ref position in read and record cigar op
        resolve_cigar(&mut r, pos);
        let qual = *r.rec.qual().get(r.qpos).unwrap_or(&0);

        if qual < min_baseq {
            drop(r);
            rbuf.backup_buf.push(raw);
            continue;
        }

        out.intake(&r, *ref_sequence)?;

        drop(r);
        rbuf.backup_buf.push(raw);
    }

    rbuf.reset();
    Ok(generated)
}
