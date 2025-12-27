use crate::{
    alignment::PileupAlignmentRef,
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

pub const UNINIT_POS: i64 = i64::MAX - 1;
pub const UNINIT_TID: i32 = i32::MAX - 1;

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
        if r.rec.pos() > pos || r.rec.tid() > tid {
            drop(r);
            rbuf.backup_buf.push(raw);
            skip_remainder_of_buf = true;
            continue;
        }

        // record is old and no longer overlaps the query coordinate. Discard.
        if read_ends_before_pos(&r, pos) {
            rbuf.depth -= 1;
            continue;
        }

        if r.rec.tid() < tid {
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

/// The iterator responsible for coordinate-wise traversal of a BAM reference/region.
/// Maintains a buffer of pileups, iterator state (current reference and position), and
/// auxiliary structs needed for variant detection and output.
pub struct PileupIterator<T: OrderedPileupOutput> {
    tid: i32,
    next_tid: i32,
    pos: i64,
    next_pos: i64,
    pub max_pos: i64,
    rbuf: ReadBuffer,
    output: Option<T>,
    dest: OutputMethod<T>,
    pub reader: BamReader,
    refseq: Option<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    show_all: bool,
    realign: bool,
    min_baseq: u8,
    redo_baq: bool,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

impl<T: OrderedPileupOutput + 'static> PileupIterator<T> {
    pub fn new(src: &BamDataSource, params: &PileupParams, output: T, dest: OutputMethod<T>) -> Result<Self, Error> {
        let reader = BamReader::new(src, 2)?;

        let rbuf = ReadBuffer::new(params.depth, params.disable_overlaps);

        let read_filter = ReadFilter::new(
            params.min_mapq,
            params.count_orphans,
            params.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        let cur_rec = Record::new();

        let tid @ next_tid = UNINIT_TID;
        let pos @ next_pos @ max_pos = UNINIT_POS;

        let show_all = params.show_empty_coords;
        let min_baseq = params.min_baseq;

        let refseq = if let Some(ref_file) = &params.refseq {
            Some(RefSeq::from_file(ref_file)?)
        } else {
            None
        };

        Ok(Self {
            tid,
            next_tid,
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
            show_all,
            realign: !params.no_baq,
            redo_baq: params.redo_baq,
        })
    }

    pub fn init_to_region(&mut self, reg: &GenomeInterval) -> Result<(), Error> {
        self.tid = i32::try_from(reg.tid)?;
        self.pos = reg.start;
        self.next_pos = reg.start;
        self.max_pos = reg.end;
        self.reader.init_to_ref(self.tid as u32, self.pos, self.max_pos)?;

        if let Some(refseq) = &mut self.refseq {
            refseq.load_seq(&self.reader.cur_ref)?;
        }
        // self.init_to_ref(false)?;

        Ok(())
    }

    /// Generate a pileup from all bases passing the minimum quality filter and covering the
    /// iterator's current reference position.
    ///
    /// If allocate is true, allocate a new output type T
    #[inline(always)]
    pub fn set_pileup(&mut self) -> Result<(), Error> {
        assert!(self.rbuf.backup_buf.is_empty());

        // don't bother going through read buffer if it starts beyond the
        // current coordinate
        if self.rbuf.start() > self.pos {
            return Ok(());
        }

        let dest = &mut self.dest;

        let ref_sequence = &mut self.refseq.as_ref().map(|r| r.yield_seq());
        let rbuf = &mut self.rbuf;

        match dest {
            OutputMethod::WriteDirectly(ref mut writer) => {
                let mut output = self.output.take().unwrap();
                output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);
                let generated = generate_pileup(rbuf, ref_sequence, &mut output, self.pos, self.tid, self.min_baseq)?;
                if generated || output.depth() > 0 || self.show_all {
                    output.write(writer)?;
                } else {
                    output.clear();
                }
                self.output = Some(output);
            }

            OutputMethod::QueueForOutput(output_chunk) => {
                let output = output_chunk.get_current_mut();
                output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);
                let generated = generate_pileup(rbuf, ref_sequence, output, self.pos, self.tid, self.min_baseq)?;
                if generated || output.depth() > 0 || self.show_all {
                    output_chunk.advance();
                } else {
                    output_chunk.tombstone();
                }
            }
        }

        Ok(())
    }

    // This is only called when we are iterating over the entire bam and we encounter a read at
    // another reference.
    pub fn increment_tid(&mut self) -> Result<(), Error> {
        self.tid = self.next_tid;
        self.pos = self.next_pos;
        self.reader.init_to_ref(self.tid as u32, 0, i64::MAX)?;

        if let Some(refseq) = &mut self.refseq {
            refseq.load_seq(&self.reader.cur_ref)?;
        }

        Ok(())
    }

    pub fn auto_loop(&mut self, region: &GenomeInterval) -> Result<(), Error> {
        self.init_to_region(region)?;

        loop {
            match self.intake()? {
                IterResult::Generated => {
                    while self.pos < self.next_pos {
                        self.set_pileup()?;
                        self.pos += 1;
                    }
                }
                IterResult::ReferenceEnd => {
                    while self.rbuf.depth > 1 {
                        self.set_pileup()?;
                        self.pos += 1;
                    }

                    self.increment_tid()?;
                }
                IterResult::NoData => break,
            }
        }

        // process what's left of the buffer after we've hit the end
        while self.rbuf.depth > 0 {
            self.set_pileup()?;
            self.pos += 1;
        }

        // if we are storing output in intermediate buffer, flush it.
        match &mut self.dest {
            OutputMethod::WriteDirectly(_) => (),
            OutputMethod::QueueForOutput(out) => out.flush(),
        }

        Ok(())
    }

    pub fn realign(&mut self, plp_ref: PileupAlignmentRef) -> Result<i32, Error> {
        if let Some(refseq) = &mut self.refseq {
            if self.realign {
                let rec = &mut plp_ref.borrow_mut().rec;
                let flag = if self.redo_baq { 7 } else { 3 };
                realign_record(rec, refseq.yield_seq(), refseq.len(), flag)?;
            }
        }

        Ok(0)
    }

    // load the read buffer until we either 1) run out of data or 2) hit a read at the next
    // position/tid.
    pub fn intake(&mut self) -> Result<IterResult, Error> {
        while self.pos == self.next_pos {
            // we need to keep reading until we have gathered all reads overlapping a position.
            //
            // TODO: move the IO reading logic outside
            if let Some(read) = self.reader.read_no_alloc(&mut self.cur_rec) {
                read?;
                let r = &self.cur_rec;

                if r.is_unmapped() {
                    continue;
                }

                if !self.read_filter.check_read(r) {
                    continue;
                }

                if r.pos() > self.max_pos {
                    break;
                }

                let ret = self.rbuf.attempt_push(r, self.pos, self.tid);

                match ret {
                    BufPushResult::Unmapped => panic!(),

                    BufPushResult::DifferentReference(pushed) => {
                        self.next_tid = r.tid();
                        self.next_pos = r.pos();
                        self.realign(pushed)?;
                        return Ok(IterResult::ReferenceEnd);
                    }

                    BufPushResult::Pushed(pushed) => {
                        self.next_pos = r.pos();
                        self.realign(pushed)?;
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
                return Ok(IterResult::NoData);
            }
        }

        // we're at the same TID as when we started, but we hit a read starting at the next
        // coordinate.
        Ok(IterResult::Generated)
    }
}
