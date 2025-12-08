use crate::{
    bamio::{BamDataSource, BamReader},
    cigar_resolve::resolve_cigar,
    output::{OrderedPileupOutput, OutputMethod},
    params::PileupParams,
    position_queue::{GenomeInterval, PositionQueue},
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
        if r.rec.pos() > pos {
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

        generated = true;

        // advance to the current ref position in read and record cigar op
        resolve_cigar(&mut r, pos);
        let qual = *r.rec.qual().get(r.qpos).unwrap_or(&0);

        if qual < min_baseq {
            drop(r);
            rbuf.backup_buf.push(raw);
            continue;
        }

        // self.output.intake(&r, ref_sequence)?;
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
pub struct PileupIterator<T: OrderedPileupOutput, W: std::io::Write> {
    tid: i32,
    pos: i64,
    next_pos: i64,
    pub max_pos: i64,
    rbuf: ReadBuffer,
    output: Option<T>,
    dest: OutputMethod<W, T>,
    pub reader: BamReader,
    refseq: Option<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    show_all: bool,
    min_baseq: u8,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

impl<T: OrderedPileupOutput + 'static, W: std::io::Write> PileupIterator<T, W> {
    pub fn new(src: &BamDataSource, params: &PileupParams, output: T, dest: OutputMethod<W, T>) -> Result<Self, Error> {
        let reader = BamReader::new(src, 2)?;

        let rbuf = ReadBuffer::new(params.depth, params.disable_overlaps);

        let read_filter = ReadFilter::new(
            params.min_mapq,
            params.count_orphans,
            params.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        let cur_rec = Record::new();

        let tid = UNINIT_TID;
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
        })
    }

    pub fn init_to_region(&mut self, reg: &GenomeInterval) -> Result<(), Error> {
        self.tid = i32::try_from(reg.tid)?;
        self.pos = reg.start;
        self.next_pos = reg.start;
        self.max_pos = reg.end;
        self.init_to_ref(false)?;

        Ok(())
    }

    // pub fn _auto_loop_step(&mut self, queue: &PositionQueue, step: usize) {}

    /// Run the iterator over the entire region without interruption.
    pub fn _auto_loop(&mut self, queue: &PositionQueue) -> Result<(), Error> {
        for reg in &queue.queue {
            self.init_to_region(reg)?;

            loop {
                match self.next()? {
                    IterResult::NoData | IterResult::ReferenceEnd => break,
                    IterResult::Generated => continue,
                }
            }

            match &mut self.dest {
                OutputMethod::WriteDirectly(_) => (),
                OutputMethod::QueueForOutput(sender, outputs) => {
                    for o in outputs.drain(..) {
                        sender.send(o)?;
                    }
                }
            }
        }

        Ok(())
    }

    pub fn _auto_loop_yield_batch(mut self, queue: &PositionQueue) -> Result<Vec<T>, Error> {
        assert_eq!(queue.len(), 1);
        self.init_to_region(&queue.queue[0])?;

        loop {
            match self.next()? {
                IterResult::NoData | IterResult::ReferenceEnd => break,
                IterResult::Generated => continue,
            }
        }

        match &mut self.dest {
            OutputMethod::WriteDirectly(_) => anyhow::bail!("Cannot output vec of reads when we output them directly"),
            OutputMethod::QueueForOutput(_sender, outputs) => {
                let out = std::mem::take(outputs);
                Ok(out)
            }
        }
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
                let generated = generate_pileup(rbuf, ref_sequence, &mut output, self.pos, self.min_baseq)?;
                if generated || output.depth() > 0 || self.show_all {
                    output.write(writer)?;
                }
                self.output = Some(output);
            }

            OutputMethod::QueueForOutput(_sender, output_chunk) => {
                let mut output = T::new();
                output.set_ref_info(self.tid, self.pos, &self.reader.cur_ref, *ref_sequence);
                let generated = generate_pileup(rbuf, ref_sequence, &mut output, self.pos, self.min_baseq)?;
                if generated || output.depth() > 0 || self.show_all {
                    output_chunk.push(output);
                }
            }
        }

        Ok(())
    }

    pub fn init_to_ref(&mut self, inc: bool) -> Result<IterResult, Error> {
        // todo: check if this works for bam files without refs in header
        if self.tid == UNINIT_TID {
            self.tid = 0;
        } else if inc {
            self.tid += 1;
        }

        self.reader.init_to_ref(self.tid as u32, self.pos, self.max_pos)?;

        if let Some(refseq) = &mut self.refseq {
            refseq.load_seq(&self.reader.cur_ref)?;
        };

        self.next_pos = 0;

        Ok(IterResult::Generated)
    }

    #[inline(always)]
    pub fn next(&mut self) -> Result<IterResult, Error> {
        while self.pos < self.next_pos {
            self.set_pileup()?;
            self.pos += 1;
        }

        while let Some(read) = self.reader.read_no_alloc(&mut self.cur_rec) {
            read?;
            let r = &self.cur_rec;

            if r.is_unmapped() {
                continue;
            }

            if r.pos() < self.pos {
                continue;
            }

            if !self.read_filter.check_read(r) {
                continue;
            }

            if r.tid() != self.tid {
                panic!();
            }

            if r.pos() > self.max_pos {
                break;
            }

            let ret = self.rbuf.attempt_push(r, self.pos, self.tid);

            match ret {
                BufPushResult::Unmapped => panic!(),
                BufPushResult::DifferentReference => return Ok(IterResult::ReferenceEnd),

                // if we hit depth limit, we exhaust all reads at this position to avoid dealing
                // with them at self.pos + 1.
                BufPushResult::MaxDepthMet => continue,
                BufPushResult::Pushed => self.next_pos = r.pos(),
            }

            if self.next_pos != self.pos && self.next_pos <= self.max_pos {
                return Ok(IterResult::Generated);
            }

            if self.next_pos > self.max_pos {
                break;
            }
        }

        while self.pos < self.max_pos {
            self.set_pileup()?;
            self.pos += 1;
        }

        Ok(IterResult::NoData)
    }
}
