use crate::{
    bamio::{BamDataSource, BamReader},
    baq::realign_record,
    cigar_resolve::resolve_cigar,
    engine::MIN_BAM_READ_THREADS,
    output::{OrderedPileupOutput, OutputMethod},
    overlap::MapOverlaps,
    params::PileupParams,
    position_queue::GenomeInterval,
    read_buf::{BufPushResult, ReadBuffer},
    read_filter::ReadFilter,
    refseq::RefSeq,
    utils::read_ends_before_pos,
};

use anyhow::Error;
use rust_htslib::bam::Record;

#[derive(Clone)]
enum EmitStrategy {
    /// output absolutely nothing
    #[allow(dead_code)]
    Nothing,

    /// output for a position if it has coverage
    ByPos,

    /// output for an entire reference if it has coverage anywhere
    ByRef,

    /// output all coordinates for all references regardless of coverage
    Everything,
}

pub struct PileupIterator<T: OrderedPileupOutput> {
    tid: i32,
    next_tid: i32,
    last_tid_with_cov: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,

    emit: EmitStrategy,

    rbuf: ReadBuffer,
    output: Option<T>,
    dest: OutputMethod<T>,
    pub reader: BamReader,
    refseq: Option<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    realign: bool,
    min_baseq: u8,
    min_mapq: u8,
    redo_baq: bool,

    read_len: usize,
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
        let reader = BamReader::new(src, MIN_BAM_READ_THREADS)?;

        let rbuf = ReadBuffer::new(params.depth, params.disable_overlaps);

        let read_filter = ReadFilter::new(
            params.count_orphans,
            params.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        let cur_rec = Record::new();

        let emit = if params.show_empty_regions {
            EmitStrategy::Everything
        } else if params.show_empty_coords {
            EmitStrategy::ByRef
        } else {
            EmitStrategy::ByPos
        };

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
            pos,
            next_pos,
            max_pos,
            rbuf,
            output: Some(output),
            dest,
            emit,
            reader,
            read_filter,
            refseq,
            cur_rec,
            min_baseq,
            min_mapq,
            realign: !params.no_baq,
            redo_baq: params.redo_baq,
            read_len: 0,
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
        if matches!(self.emit, EmitStrategy::ByPos) {
            if self.rbuf.head.tid == self.tid && self.rbuf.head.pos > self.pos {
                skip = true;
            }

            if self.rbuf.head.tid > self.tid {
                skip = true
            }
        }

        let ref_sequence = &mut self.refseq.as_ref().and_then(|r| r.yield_seq());
        let rbuf = &mut self.rbuf;
        let generated;
        let depth;

        {
            let output = self.dest.cur();
            output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);

            if !skip {
                generated = generate_pileup(rbuf, ref_sequence, output, self.pos, self.tid, self.min_baseq)?;
                depth = output.depth();
            } else {
                generated = false;
                depth = 0;
            };
        }

        let written = match self.emit {
            EmitStrategy::Nothing => self.dest.reject()?,
            EmitStrategy::ByPos => self.dest.check(generated || depth > 0)?,
            EmitStrategy::ByRef => self
                .dest
                .check((self.tid == self.last_tid_with_cov) || self.rbuf.head.tid == self.tid)?,
            EmitStrategy::Everything => self.dest.take()?,
        };

        if written {
            self.last_tid_with_cov = self.tid;
        }

        Ok(())
    }

    fn set_ref(&mut self, interval: GenomeInterval) -> Result<(), Error> {
        if interval.tid >= self.reader.header.target_count() as i64 {
            anyhow::bail!("Interval has TID exceeding header maximum!");
        }

        let mut output = self.output.take().unwrap();

        // purge read buffer to remove any reads spanning the old ref to update head and tail.
        generate_pileup(
            &mut self.rbuf,
            &self.refseq.as_ref().and_then(|r| r.yield_seq()),
            &mut output,
            i64::MAX,
            self.tid,
            self.min_baseq,
        )?;

        output.clear();
        self.output = Some(output);

        self.reader
            .init_to_ref(interval.tid as u32, interval.start, interval.end)?;
        self.pos = interval.start;
        self.next_pos = interval.start;

        self.tid = interval.tid as i32;
        self.next_tid = self.tid;

        self.max_pos = interval.end - 1;

        if let Some(refseq) = &mut self.refseq {
            refseq.load_seq(&self.reader.cur_ref)?;
        }

        Ok(())
    }

    // load the read buffer until we either 1) run out of data or 2) hit a read at the next
    // position/tid.
    #[inline(always)]
    pub fn intake(&mut self) -> Result<IterResult, Error> {
        if self.reader.eof {
            return Ok(IterResult::ReferenceEnd);
        }

        loop {
            // we need to keep reading until we have gathered all reads overlapping a position.
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
                    self.rbuf.attempt_push(self.tid, self.pos, r)?;
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

                        if self.next_pos > self.pos {
                            return Ok(IterResult::Generated);
                        }
                    }

                    // if we've capped our buffer to a given depth, we'll iterate over all
                    // remaining reads spanning this coordinate before stopping to generate
                    // pileups. This way we won't have to deal with them at the next position.
                    BufPushResult::MaxDepthMet | BufPushResult::BeforePos => {
                        continue;
                    }
                }
            } else {
                // we ran out of reads.
                self.reader.eof = true;
                return Ok(IterResult::ReferenceEnd);
            }
        }
    }

    pub fn auto_loop2(&mut self, intervals: &[GenomeInterval]) -> Result<(), Error> {
        self.read_len = BamReader::sample_read_len(&self.reader.src)?;

        for interval in intervals {
            self.set_ref(interval.clone())?;
            self.process_single_ref()?;

            match self.emit {
                EmitStrategy::ByRef | EmitStrategy::Everything => {
                    while self.pos <= self.max_pos {
                        self.set_pileup()?;
                        self.pos += 1;
                    }
                }
                _ => (),
            }
        }
        Ok(())
    }

    pub fn process_single_ref(&mut self) -> Result<(), Error> {
        loop {
            // eprintln!(
            //     "MAIN: {} {} {} {}",
            //     self.pos, self.next_pos, self.max_pos, self.rbuf.depth
            // );
            match self.intake()? {
                IterResult::Generated => {
                    // eprintln!("GENERATED: {} {} {}", self.pos, self.rbuf.head.tid, self.rbuf.head.pos);
                    if self.rbuf.head.pos < self.pos
                        || (matches!(self.emit, EmitStrategy::ByRef) && self.rbuf.head.tid == self.tid)
                        || matches!(self.emit, EmitStrategy::Everything)
                    {
                        while self.pos < self.next_pos && self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1;
                        }
                    } else {
                        self.pos = self.rbuf.head.pos;

                        while self.pos < self.next_pos && self.pos <= self.max_pos {
                            self.set_pileup()?;
                            self.pos += 1
                        }
                    }
                }

                IterResult::ReferenceEnd => {
                    // eprintln!(
                    //     "REF END:{} {} {} {}",
                    //     self.pos, self.rbuf.head.tid, self.rbuf.head.pos, self.rbuf.depth
                    // );
                    // if we have reads for current ref still in buffer, process them until they no
                    // longer overlap with cur pos.
                    while self.rbuf.head.tid == self.tid && self.pos <= self.max_pos {
                        self.set_pileup()?;
                        self.pos += 1;
                    }

                    break;
                }
            }
        }

        // if we are storing output in intermediate buffer, flush it.
        match &mut self.dest {
            OutputMethod::WriteDirectly(_, _) => (),
            OutputMethod::QueueForOutput(out) => out.write_all()?,
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
            if let Some(ref mut overlap) = rbuf.overlap_map {
                overlap.delete_read(&r.rec);
            }
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
