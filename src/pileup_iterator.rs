use crate::{
    bamio::{BamDataSource, BamReader},
    baq::realign_record,
    cigar_resolve::resolve_cigar,
    engine::MIN_BAM_READ_THREADS,
    output::{OrderedPileupOutput, OutputFormat},
    overlap::MapOverlaps,
    params::PileupParams,
    position_queue::GenomeInterval,
    read_buf::{BufPushResult, ReadBuffer},
    read_filter::ReadFilter,
    refseq::RefSeqHandle,
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

/// The state machine performing the pileup generation and advancing across region coordinates.
/// The type T dictates what output it generates.
pub struct PileupIterator<T: OrderedPileupOutput> {
    tid: i32,
    next_tid: i32,
    last_tid_with_cov: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,

    emit: EmitStrategy,

    rbuf: ReadBuffer,
    dest: OutputFormat<T>,
    pub reader: BamReader,
    refseq: RefSeqHandle,
    read_filter: ReadFilter,
    cur_rec: Record,
    realign: bool,
    min_baseq: u8,
    min_mapq: u8,
    redo_baq: bool,

    read_len: usize,
}

impl<T: OrderedPileupOutput> PileupIterator<T> {
    /// Create a new pileup iterator from a data source (e.g. bam file), a set of query regions,
    /// input params and an output type.
    pub fn new(
        src: &BamDataSource,
        refseq: RefSeqHandle,
        params: &PileupParams,
        dest: OutputFormat<T>,
    ) -> Result<Self, Error> {
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
    /// iterator's current reference position. Importantly, generate_pileup() is where stale reads no longer
    /// overlapping the query position are removed.
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

        // let ref_sequence = &mut self.refseq.as_ref().and_then(|r| r.yield_seq());
        let rbuf = &mut self.rbuf;
        let generated;
        let depth;

        {
            let output = self.dest.cur();
            output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, &self.refseq);

            if !skip {
                generated = generate_pileup(rbuf, &self.refseq, output, self.pos, self.tid, self.min_baseq)?;
                depth = output.depth();
            } else {
                generated = false;
                depth = 0;
            };
        }

        let written = match self.emit {
            EmitStrategy::Nothing => self.dest.reject(),
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

    /// When given a region not starting at zero, rewind by 2X read length in order to populate the
    /// overlap map to ensure read mates get nulled. This should only be used in a multi-threading
    /// context.
    fn preload_region(&mut self, interval: &GenomeInterval) -> Result<(), Error> {
        let rewind = (interval.start - (2 * self.read_len) as i64).max(0);

        self.pos = rewind;
        self.next_pos = self.pos;
        self.max_pos = interval.start - 1;

        self.tid = interval.tid as i32;
        self.next_tid = self.tid;

        self.reader.init_to_ref(interval.tid as u32, self.pos, interval.end)?;

        let preset = self.emit.clone();
        self.emit = EmitStrategy::Nothing;
        self.process_single_ref()?;

        assert!(self.pos <= interval.start);
        self.pos = interval.start;
        self.next_pos = self.pos;

        self.emit = preset;

        Ok(())
    }

    /// Update iterator state and prepare ref-specific data given a new interval.
    fn set_ref(&mut self, interval: GenomeInterval) -> Result<(), Error> {
        if interval.tid >= self.reader.header.target_count() as i64 {
            anyhow::bail!("Interval has TID exceeding header maximum!");
        }

        let output = self.dest.cur();

        // purge read buffer to remove any reads spanning the old ref to update head and tail.
        generate_pileup(&mut self.rbuf, &self.refseq, output, i64::MAX, self.tid, self.min_baseq)?;

        output.clear();

        if interval.start != 0 && self.rbuf.overlap_map.is_some() {
            self.preload_region(&interval)?;
        } else {
            self.reader
                .init_to_ref(interval.tid as u32, interval.start, interval.end)?;
            self.pos = interval.start;
            self.next_pos = interval.start;
        }

        self.tid = interval.tid as i32;
        self.next_tid = self.tid;

        self.max_pos = interval.end - 1;

        // if let Some(refseq) = &mut self.refseq {
        //     refseq.load_seq(&self.reader.cur_ref)?;
        // }

        Ok(())
    }

    /// Read from an input BAM until we find either 1) read starting at a coordinate beyond the end of the
    /// read buffer, 2) run out of data, 3) exceed the max position of this region, or 4) encounter
    /// a read mapping to a new reference.
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

                if self.realign {
                    if let Some(refseq) = self.refseq.as_ref() {
                        let flag = if self.redo_baq { 7 } else { 3 };
                        realign_record(r, refseq, refseq.len() as i64, flag)?;
                    }
                }

                // we passed queried region
                if r.pos() > self.max_pos {
                    self.rbuf.attempt_push(self.tid, self.pos, r)?;
                    return Ok(IterResult::ReferenceEnd);
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

    /// main function of the PileupIterator: run it on all the query intervals given.
    pub fn auto_loop2(&mut self, interval: &GenomeInterval) -> Result<(), Error> {
        self.read_len = BamReader::sample_read_len(&self.reader.src)?;

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

        Ok(())
    }

    /// Iterate across the coordinates of a query interval and produce output.
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

        Ok(())
    }
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
}

#[inline(always)]
/// Perform the pileup given a read buffer, optional ref sequence, an output destination, and query
/// position. Importantly, reads found to no longer overlap (pos, tid) will be removed from the
/// buffer.
pub fn generate_pileup<T: OrderedPileupOutput>(
    rbuf: &mut ReadBuffer,
    refseq: &RefSeqHandle,
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
            drop(r);
            rbuf.backup_buf.push(raw);
            skip_remainder_of_buf = true;
            continue;
        }

        // record is old and no longer overlaps the query coordinate. We discard it by not adding
        // it to the alternate buffer.
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

        out.intake(&r, refseq)?;

        drop(r);
        rbuf.backup_buf.push(raw);
    }

    rbuf.reset();
    Ok(generated)
}
