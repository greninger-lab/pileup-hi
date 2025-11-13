use crate::{
    bamio::{BamDataSource, BamReader},
    cigar_resolve::resolve_cigar,
    params::PileupParams,
    pileup_string::PileupString,
    position_queue::PositionQueue,
    read_buf::{BufPushResult, ReadBuffer},
    read_filter::ReadFilter,
    refseq::RefSeq,
    utils::read_ends_before_pos,
};

use crossbeam::channel::Sender;

use anyhow::Error;
use rust_htslib::bam::Record;

pub const UNINIT_POS: i64 = i64::MAX - 1;
pub const UNINIT_TID: i32 = i32::MAX - 1;

/// The iterator responsible for coordinate-wise traversal of a BAM reference/region.
/// Maintains a buffer of pileups, iterator state (current reference and position), and
/// auxiliary structs needed for variant detection and output.
pub struct PileupIterator {
    tid: i32,
    pos: i64,
    next_pos: i64,
    pub max_pos: i64,
    rbuf: ReadBuffer,
    pileup_writer: PileupString,
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

impl PileupIterator {
    pub fn new(
        src: &BamDataSource,
        params: &PileupParams,
        out_handle: Option<Sender<PileupString>>,
    ) -> Result<Self, Error> {
        let reader = BamReader::new(src, num_cpus::get())?;

        let pileup_writer = PileupString::new();
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
            pileup_writer,
            reader,
            read_filter,
            refseq,
            cur_rec,
            min_baseq,
            show_all,
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

    /// Generate a pileup from all bases passing the minimum quality filter and covering the
    /// iterator's current reference position.
    pub fn set_pileup(&mut self) -> Result<bool, Error> {
        assert!(self.rbuf.backup_buf.is_empty());

        let (ref_sequence, ref_base) = if let Some(refseq) = &self.refseq {
            let r = refseq.yield_seq();
            (Some(r), *r.get(self.pos as usize).unwrap_or(&b'N'))
        } else {
            (None, b'N')
        };

        // don't bother going through read buffer if it starts beyond the
        // current coordinate
        if self.rbuf.start() > self.pos {
            return Ok(false);
        }

        let mut generated = false;
        let mut skip_remainder_of_buf = false;

        self.pileup_writer
            .update(self.tid, self.pos, ref_base, &self.reader.cur_ref);

        for raw in self.rbuf.rbuf.drain(..) {
            // from a previous record, we decided to skip all remaining records in this buffer.
            if skip_remainder_of_buf {
                self.rbuf.backup_buf.push(raw);
                continue;
            }

            let mut r = raw.borrow_mut();

            // record starts beyond position, which means that the remainder of the buffer does
            // too. Skip the rest of the records.
            if r.rec.pos() > self.pos {
                drop(r);
                self.rbuf.backup_buf.push(raw);
                skip_remainder_of_buf = true;
                continue;
            }

            // record is old and no longer overlaps the query coordinate. Discard.
            if read_ends_before_pos(&r, self.pos) {
                self.rbuf.depth -= 1;
                continue;
            }

            generated = true;

            // advance to the current ref position in read and record cigar op
            resolve_cigar(&mut r, self.pos);
            let qual = *r.rec.qual().get(r.qpos).unwrap_or(&0);

            if qual < self.min_baseq {
                drop(r);
                self.rbuf.backup_buf.push(raw);
                continue;
            }

            self.pileup_writer.intake(&r, ref_sequence)?;

            drop(r);
            self.rbuf.backup_buf.push(raw);
        }

        if self.pileup_writer.depth > 0 || self.show_all {
            self.pileup_writer.write()?;
        }

        self.rbuf.reset();

        Ok(generated)
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

//#[cfg(test)]
//mod tests {
//    use super::*;
//    use crate::alignment::CigarState;
//    use rust_htslib::bam::{record::CigarString, Record};

//    #[test]
//    pub fn cig_test1() {
//        let cig = Vec::from([Cigar::Match(76)]);
//        assert_eq!(cig[0].len(), 76)
//    }

//    #[test]
//    pub fn cig_test2() {
//        let mut record = Record::new();
//        record.set(
//            b"read1",
//            Some(&CigarString(vec![Cigar::Match(4), Cigar::Equal(1)])),
//            b"AAAAG",
//            b"#####",
//        );

//        record.set_pos(1);

//        let mut cstate = CigarState {
//            cig: record.cigar(),
//            icig: 0,
//            iseq: 0,
//            bam_pos: 1,
//            qpos: 0,
//            del: false,
//            read_len_from_cigar: 0,
//        };

//        let mut ret = cigar_get_pos(&mut cstate, 4);
//        assert_eq!(ret, Some(Cigar::Match(4)));
//        ret = cigar_get_pos(&mut cstate, 5);
//        assert_eq!(ret, Some(Cigar::Match(1)))
//    }

//    #[test]
//    pub fn cig_test3() {
//        let mut record = Record::new();
//        record.set(
//            b"read1",
//            Some(&CigarString(vec![
//                Cigar::Match(4),
//                Cigar::Equal(1),
//                Cigar::Ins(2),
//                Cigar::Match(3),
//            ])),
//            b"AAAAGTTTTT",
//            b"##########",
//        );

//        record.set_pos(104);

//        let mut cstate = CigarState {
//            cig: record.cigar(),
//            icig: 0,
//            iseq: 0,
//            bam_pos: 104,
//            qpos: 0,
//            del: false,
//            read_len_from_cigar: 0,
//        };

//        let mut ret = cigar_get_pos(&mut cstate, 107);
//        assert_eq!(ret, Some(Cigar::Match(4)));

//        ret = cigar_get_pos(&mut cstate, 108);
//        assert_eq!(ret, Some(Cigar::Ins(2)));

//        ret = cigar_get_pos(&mut cstate, 109);
//        assert_eq!(ret, Some(Cigar::Match(3)));
//    }

//    #[test]
//    pub fn cig_test4() {
//        let mut record = Record::new();
//        record.set(
//            b"read1",
//            Some(&CigarString(vec![Cigar::Match(1), Cigar::Del(4), Cigar::Match(3)])),
//            b"AATTTT",
//            b"##EEEE",
//            //012345
//        );

//        record.set_pos(104);

//        let mut cstate = CigarState {
//            cig: record.cigar(),
//            icig: 0,
//            iseq: 0,
//            bam_pos: 104,
//            qpos: 0,
//            del: false,
//            read_len_from_cigar: 0,
//        };

//        let mut ret = cigar_get_pos(&mut cstate, 104);
//        assert_eq!(ret, Some(Cigar::Del(4)));

//        ret = cigar_get_pos(&mut cstate, 105);
//        assert_eq!(ret, Some(Cigar::Del(4)));
//        assert_eq!(cstate.qpos, 1);

//        ret = cigar_get_pos(&mut cstate, 106);
//        assert_eq!(ret, Some(Cigar::Del(4)));
//        assert_eq!(cstate.qpos, 1);

//        ret = cigar_get_pos(&mut cstate, 107);
//        assert_eq!(ret, Some(Cigar::Del(4)));
//        assert_eq!(cstate.qpos, 1);

//        ret = cigar_get_pos(&mut cstate, 108);
//        assert_eq!(ret, Some(Cigar::Del(4)));
//        assert_eq!(cstate.qpos, 1);

//        ret = cigar_get_pos(&mut cstate, 109);
//        assert_eq!(ret, Some(Cigar::Match(3)));
//        assert_eq!(cstate.qpos, 1);

//        ret = cigar_get_pos(&mut cstate, 110);
//        assert_eq!(ret, Some(Cigar::Match(3)));
//        assert_eq!(cstate.qpos, 2);
//    }
//}
