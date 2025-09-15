use crate::bamio::{BamReader, BamWriter};
use crate::params::Params;
use crate::pileup_writer::PileupWriter;
use crate::position_queue::PositionQueue;
use crate::read_buf;
use crate::read_filter::ReadFilter;
use crate::realigner::{AlignerReference, Realigner};
use crate::refseq::RefSeq;
use crate::utils::{cigar_get_pos, read_ends_before_pos};

use anyhow::{Context, Error};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use std::cell::RefCell;

pub const UNINIT_POS: i64 = i64::MAX - 1;
pub const UNINIT_TID: i32 = i32::MAX - 1;

/// A counter class to track the number of reads differing from the reference at a given pileup.
/// Holds the qualities of all bases, and the number of reads with ref-matching vs differing bases.
///
/// [nmatch] = number of reads consuming reference at this positiion, which includes nalt.
/// [nalt] = number of nmatch reads that differ from reference sequence
/// [ndel] = number of reads with deletion following this position
/// [nins] = number of reads with insertion following this position
pub struct PileupPosition {
    quals: Vec<u8>,
    nmatch: u32,
    nalt: u32,
    ndel: u32,
    nins: u32,
    ref_base: u8,
    register: fn(&mut Self, u8, u8),
}

impl PileupPosition {
    fn _uptake(&mut self, qual: u8, base: u8) {
        self.quals.push(qual);

        match base == self.ref_base {
            true => self.nmatch += 1,
            false => self.nalt += 1,
        }
    }

    fn _ignore(&mut self, _qual: u8, base: u8) {
        match base == self.ref_base {
            true => self.nmatch += 1,
            false => self.nalt += 1,
        }
    }

    pub fn new(gather_align_metrics: bool) -> Self {
        let uptake_method = match gather_align_metrics {
            true => PileupPosition::_uptake,
            false => PileupPosition::_ignore,
        };

        Self {
            quals: Vec::new(),
            nmatch: 0,
            nalt: 0,
            ndel: 0,
            nins: 0,
            ref_base: b'N',
            register: uptake_method,
        }
    }

    /// Is a current pileup position in need of reassembly, i.e. is there a significant amount of
    /// differences (indels, substitutions) relative to reference?
    pub fn is_active(&mut self, floor: f32, ceil: f32, denom: f32) -> Option<f32> {
        let ratio = (self.nalt + self.ndel + self.nins) as f32 / denom;
        match (ratio >= floor && ratio <= ceil) && (self.nins + self.ndel) as f32 / denom > 0.1 {
            true => Some(ratio),
            false => None,
        }
    }

    pub fn depth(&mut self) -> u32 {
        self.ndel + self.nins + self.nmatch + self.nalt
    }

    pub fn clear(&mut self) {
        self.quals.clear();
        self.nalt = 0;
        self.nins = 0;
        self.ndel = 0;
        self.nmatch = 0;
    }
}

/// The iterator responsible for coordinate-wise traversal of a BAM reference/region.
/// Maintains a buffer of pileups, iterator state (current reference and position), and
/// auxiliary structs needed for variant detection and output.
pub struct PileupIterator {
    tid: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,
    store: PileupPosition,
    tid_count: i32,
    rbuf: read_buf::ReadBuffer,
    realigner: Option<Realigner>,
    show_all: bool,
    pileup_writer: PileupWriter,
    pub reader: BamReader,
    refseq: RefCell<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    min_baseq: u8,
    read_discard: RefCell<BamWriter>,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

impl PileupIterator {
    pub fn new(params: &Params) -> Result<Self, Error> {
        let tid = params.inp.tid.unwrap_or(UNINIT_TID);
        let pos @ next_pos @ max_pos = params.inp.pos.unwrap_or(UNINIT_POS);
        let reader = BamReader::new(&params.inp)?;
        let pileup_writer = PileupWriter::new();

        let (store, realigner) = match params.plp.indel_realign {
            true => (PileupPosition::new(true), Some(Realigner::build_empty()?)),
            false => (PileupPosition::new(false), None),
        };

        let read_discard = match &params.outp.output_realigned {
            Some(output) => RefCell::new(BamWriter::new_from_template(&reader.header, output)?),
            None => RefCell::new(BamWriter::void(&reader.header)?),
        };

        let rbuf = read_buf::ReadBuffer::new(params.inp.depth, params.plp.disable_overlaps);

        let show_all = params.plp.show_empty_coords;
        let cur_rec = Record::new();
        let min_baseq = params.plp.min_baseq;
        let max_tid = reader.max_tid;

        let read_filter = ReadFilter::new(
            params.plp.min_mapq,
            params.plp.count_orphans,
            params.plp.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.plp.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        let refseq = match &params.inp.refseq {
            Some(ref_file) => RefCell::new(RefSeq::from_file(ref_file.clone())?),
            None => RefCell::new(RefSeq::new_empty()),
        };

        Ok(Self {
            tid,
            pos,
            next_pos,
            realigner,
            max_pos,
            tid_count: max_tid,
            rbuf,
            pileup_writer,
            reader,
            store,
            min_baseq,
            read_discard,
            read_filter,
            show_all,
            refseq,
            cur_rec,
        })
    }

    pub fn _auto_loop(&mut self, queue: &PositionQueue) -> Result<(), Error> {
        for reg in &queue.queue {
            self.tid = i32::try_from(reg.tid)?;
            self.pos = reg.start;
            self.next_pos = reg.start;
            self.max_pos = reg.end;
            self.init_to_ref(false)?;

            loop {
                match self.next()? {
                    IterResult::NoData | IterResult::ReferenceEnd => break,
                    IterResult::Generated => continue,
                }
            }
        }

        Ok(())
    }

    pub fn establish_position_context(&mut self) -> Result<(), Error> {
        self.store.clear();
        self.store.ref_base = self.refseq.borrow().get_base(self.pos as u64)?;

        Ok(())
    }

    /// Generate a pileup from all bases passing the minimum quality filter and covering the
    /// iterator's current reference position.
    pub fn set_pileup(&mut self) -> Result<bool, Error> {
        assert!(self.rbuf.backup_buf.is_empty());
        let mut generated = false;

        self.establish_position_context()?;
        let mut refseq = self.refseq.borrow_mut();
        let mut discard = self.read_discard.borrow_mut();

        for raw in self.rbuf.rbuf.drain(..) {
            let mut r = raw.borrow_mut();

            if read_ends_before_pos(&r.rec, self.pos) {
                discard.write_record(&r.rec)?;
                self.rbuf.depth -= 1;
                continue;
            }

            // advance to the current ref position in read and record cigar op
            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32);
            let qual = *r.rec.qual().get(r.cstate.qpos).unwrap_or(&0);
            let readbase = r.rec.seq()[r.cstate.qpos];

            if qual < self.min_baseq {
                drop(r);
                self.rbuf.backup_buf.push(raw);
                continue;
            }

            if let Some(allele) = ret {
                match allele {
                    Cigar::Match(_) => {
                        self.pileup_writer.write_match(
                            &r.rec,
                            r.cstate.qpos as u32,
                            self.pos,
                            self.store.ref_base,
                        )?;

                        (self.store.register)(&mut self.store, qual, readbase);
                    }

                    Cigar::Ins(_) => {
                        self.pileup_writer.write_match(
                            &r.rec,
                            r.cstate.qpos as u32,
                            self.pos,
                            self.store.ref_base,
                        )?;

                        self.pileup_writer.write_insertion(
                            &r.cstate,
                            &r.rec,
                            r.cstate.qpos as u32,
                        )?;

                        self.store.nins += 1;
                    }

                    Cigar::Del(l) => {
                        if !r.cstate.del {
                            self.pileup_writer.write_match(
                                &r.rec,
                                r.cstate.qpos as u32,
                                self.pos,
                                self.store.ref_base,
                            )?;

                            let del_seq = refseq.get_interval(
                                self.pos as u64 + 1,
                                (self.pos + (l as i64)) as u64,
                            )?;

                            self.pileup_writer
                                .write_deletion_start(&r.rec, del_seq, l as i64)?
                        } else {
                            self.pileup_writer.write_deletion(qual);
                        }

                        self.store.ndel += 1;
                    }

                    _ => panic!("Invalid pileup type found!"),
                }
            }

            drop(r);
            self.rbuf.backup_buf.push(raw);
        }

        let depth = self.store.depth();

        if self.show_all || depth > 0 {
            self.pileup_writer.write_pileup_str(
                self.store.ref_base,
                self.pos,
                depth,
                &self.reader.cur_ref,
            )?;
            generated = true;
        }

        if let Some(realigner) = &mut self.realigner {
            if let Some(_ratio) = self.store.is_active(0.1, 0.9, depth as f32) {
                realigner.realign_region_plp(&mut self.rbuf.backup_buf)?;
            }
        }

        self.rbuf.reset();

        Ok(generated)
    }

    pub fn init_to_ref(&mut self, inc: bool) -> Result<IterResult, Error> {
        // todo: check if this works for bam files without refs in header
        //
        if self.tid == UNINIT_TID {
            self.tid = 0;
        } else if inc {
            self.tid += 1;
        }

        self.reader.init_to_ref(self.tid as u32)?;

        let tlen = self
            .reader
            .header
            .target_len(self.tid as u32)
            .with_context(|| format!("Failed to get target length for {}", self.reader.cur_ref))?;

        let mut refseq = self.refseq.try_borrow_mut()?;

        if !refseq.is_empty() {
            // right now we just get the entire reference sequence.
            // Next step will be to load it in windows.
            refseq.load_seq(&self.reader.cur_ref, 0, tlen)?
        }

        if let Some(realigner) = &mut self.realigner {
            assert!(!refseq.is_empty());
            realigner.init_to_ref(
                AlignerReference::Sequence(refseq.yield_seq_slice()),
                Some(&self.reader.cur_ref),
            )?;
        }

        self.next_pos = 0;

        Ok(IterResult::Generated)
    }

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

            if !self.read_filter.check_read(r) {
                continue;
            }

            if r.tid() != self.tid {
                panic!();
            }

            let ret = self.rbuf.attempt_push(r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::Unmapped => panic!(),
                read_buf::BufPushResult::DifferentReference => return Ok(IterResult::ReferenceEnd),

                // if we hit depth limit, we exhaust all reads at this position to avoid dealing
                // with them at self.pos + 1.
                read_buf::BufPushResult::MaxDepthMet => continue,
                read_buf::BufPushResult::Pushed => self.next_pos = r.pos(),
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

        match self.tid + 1 == self.tid_count {
            true => Ok(IterResult::NoData),
            false => Ok(IterResult::ReferenceEnd),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarState;
    use rust_htslib::bam::{record::CigarString, Record};

    #[test]
    pub fn cig_test1() {
        let cig = Vec::from([Cigar::Match(76)]);
        assert_eq!(cig[0].len(), 76)
    }

    #[test]
    pub fn cig_test2() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(4), Cigar::Equal(1)])),
            b"AAAAG",
            b"#####",
        );

        record.set_pos(1);

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 1,
            qpos: 0,
            del: false,
        };

        let mut ret = cigar_get_pos(&mut cstate, 4);
        assert_eq!(ret, Some(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, 5);
        assert_eq!(ret, Some(Cigar::Match(1)))
    }

    #[test]
    pub fn cig_test3() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Equal(1),
                Cigar::Ins(2),
                Cigar::Match(3),
            ])),
            b"AAAAGTTTTT",
            b"##########",
        );

        record.set_pos(104);

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 104,
            qpos: 0,
            del: false,
        };

        let mut ret = cigar_get_pos(&mut cstate, 107);
        assert_eq!(ret, Some(Cigar::Match(4)));

        ret = cigar_get_pos(&mut cstate, 108);
        assert_eq!(ret, Some(Cigar::Ins(2)));

        ret = cigar_get_pos(&mut cstate, 109);
        assert_eq!(ret, Some(Cigar::Match(3)));
    }

    #[test]
    pub fn cig_test4() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(1),
                Cigar::Del(4),
                Cigar::Match(3),
            ])),
            b"AATTTT",
            b"##EEEE",
            //012345
        );

        record.set_pos(104);

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 104,
            qpos: 0,
            del: false,
        };

        let mut ret = cigar_get_pos(&mut cstate, 104);
        assert_eq!(ret, Some(Cigar::Del(4)));

        ret = cigar_get_pos(&mut cstate, 105);
        assert_eq!(ret, Some(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 106);
        assert_eq!(ret, Some(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 107);
        assert_eq!(ret, Some(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 108);
        assert_eq!(ret, Some(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 109);
        assert_eq!(ret, Some(Cigar::Match(3)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 110);
        assert_eq!(ret, Some(Cigar::Match(3)));
        assert_eq!(cstate.qpos, 2);
    }
}
