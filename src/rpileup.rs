use crate::read_buf;
use crate::read_buf::CigarState;
use anyhow::{Context, Error};
use num_cpus;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{ext::BamRecordExtensions, HeaderView, IndexedReader, Read, Record};
use std::collections::VecDeque;
use std::io::Write;

const UNINIT_POS: usize = usize::MAX - 1;
const UNINIT_TID: u32 = u32::MAX - 1;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

pub struct PileupIterator {
    tid: u32,
    pos: usize,
    next_pos: usize,
    max_pos: usize,
    show_all: bool,
    rbuf: read_buf::ReadBuffer,
    reader: IndexedReader,
    header: HeaderView,
    ref_seq: Option<Vec<u8>>,
    next_record: Option<Record>,
    coverage: u32,
    seq_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    remove_buf: VecDeque<usize>,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Pileup {
    BeforePos(),
    Op(Cigar),
    BaseEmpty(),
}

pub fn get_base_pileup(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    ref_base: u8,
) {
    let ipos = ipos as usize;
    let bam_pos = cs.bam_pos as usize;

    if pos == r.reference_end() as usize - 1 {
        seq_buf.push(LAST_POS);
        seq_buf.push(r.mapq() + 33);
    }

    if pos == bam_pos {
        seq_buf.push(FIRST_POS);
        seq_buf.push(r.mapq() + 33);
    }

    let cur_base = r.seq()[ipos];

    let mut base: u8;

    if ref_base != cur_base {
        base = cur_base;
        if r.is_reverse() {
            base.make_ascii_lowercase();
        }
    } else {
        if r.is_reverse() {
            base = R_MATCH;
        } else {
            base = F_MATCH;
        }
    };

    let cur_qual = r.qual()[ipos] + 33;

    seq_buf.push(base);
    qual_buf.push(cur_qual);
}

pub fn write_del(pos: usize, seq_buf: &mut Vec<u8>, del_len: usize) -> Result<(), Error> {
    write!(seq_buf, "-{}", del_len)?;
    for _ in pos..pos + del_len {
        seq_buf.push(b'N')
    }
    Ok(())
}

pub fn cigar_get_pos(cs: &mut CigarState, pos: u32, ipos: &mut i32) -> Pileup {
    let cig = &cs.cig;
    let ncig = cig.len();
    let mut op: Cigar;
    while cs.bam_pos <= pos {
        if cs.icig >= ncig {
            // this should never happen, since we check cigars beforehand to at least end
            // at the queried coordinate, if not pass over it.
            return Pileup::BeforePos();
        }

        op = cig[cs.icig];
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end_pos = cs.bam_pos + len - 1;

                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.iseq += len;
                    cs.icig += 1;
                    continue;
                }

                *ipos = pos as i32 - cs.bam_pos as i32 + cs.iseq as i32;
                if end_pos == pos && cs.icig + 1 < ncig {
                    let next_op = cig[cs.icig + 1];

                    match next_op {
                        Cigar::Ins(_) => return Pileup::Op(next_op),
                        Cigar::Del(_) => return Pileup::Op(next_op),
                        _ => (),
                    }
                }
                return Pileup::Op(Cigar::Match(len));
            }

            Cigar::Ins(len) | Cigar::SoftClip(len) => {
                cs.iseq += len;
                cs.icig += 1;
                continue;
            }

            Cigar::Del(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                // this coordinate comes after we already indicated the deletion, so
                // mark ipos to avoid repeating the deletion in this and subsequent plp cols
                *ipos = -1;
                return Pileup::Op(op);
            }

            Cigar::RefSkip(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                return Pileup::BaseEmpty();
            }
            _ => (),
        }
    }

    Pileup::BaseEmpty()
}

impl PileupIterator {
    pub fn new(
        bam_fname: &str,
        show_all: bool,
        tid: Option<u32>,
        pos: Option<usize>,
    ) -> Result<Self, Error> {
        let tid = tid.unwrap_or(UNINIT_TID);
        let pos @ next_pos @ max_pos = pos.unwrap_or(UNINIT_POS);
        let mut reader = IndexedReader::from_path(bam_fname)?;
        reader.set_threads(num_cpus::get())?;
        let rbuf = read_buf::ReadBuffer::new();
        let header = reader.header().clone();
        let ref_seq = None;
        let next_record = None;
        let coverage = 0;
        let remove_buf = VecDeque::with_capacity(500);
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));

        Ok(Self {
            tid,
            pos,
            next_pos,
            max_pos,
            rbuf,
            reader,
            header,
            ref_seq,
            show_all,
            next_record,
            coverage,
            seq_buf,
            qual_buf,
            remove_buf,
        })
    }

    /// Read records in to fill a read buffer spanning the current coordinate window.
    /// This will loop over records until A) a record is found that starts outside the current
    /// window, (e.g. faraway coord or different reference).
    ///
    /// When a read outside the current window is found, the [PileupIterator] will skip buffer
    /// filling / pileup generation for all coordinates between the current and the next read's
    /// start position.
    pub fn fill_buffer(&mut self) -> Result<(), Error> {
        let mut ret: read_buf::BufPushResult;
        let mut scanned = 0;

        let mut prev_pos = i64::MIN;

        if let Some(next_record) = self.next_record.take() {
            scanned += 1;
            prev_pos = next_record.pos();
            let ret = self.rbuf.push(next_record, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::AfterWindow((r, next_pos))
                | read_buf::BufPushResult::DifferentReference((r, next_pos)) => {
                    self.next_pos = self.pos + next_pos;
                    self.next_record = Some(r);
                    return Ok(());
                }
                _ => self.next_record = None,
            }
        }

        for rec in self.reader.records() {
            let r = rec?;
            scanned += 1;

            if r.is_unmapped() {
                continue;
            }

            if r.pos() < prev_pos {
                panic!("UNSORTED BAM! {} {}", r.pos(), prev_pos)
            }

            prev_pos = r.pos();

            ret = self.rbuf.push(r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::Unmapped => continue,

                read_buf::BufPushResult::AfterWindow((next_rec, next_pos))
                | read_buf::BufPushResult::DifferentReference((next_rec, next_pos)) => {
                    self.next_pos = self.pos + next_pos;
                    self.next_record = Some(next_rec);
                    break;
                }
                read_buf::BufPushResult::Pushed => scanned += 1,
            }
        }

        // if we have no reads at all to set next pos, assume
        // we've hit the end of the reference, and set next pos to MAX
        if scanned == 0 {
            self.next_pos = usize::MAX;
        }

        Ok(())
    }

    pub fn write_pileup_str(
        &mut self,
        ref_base: u8,
        nbases: usize,
        nins: usize,
        ndel: usize,
    ) -> Result<(), Error> {
        print! {"{}\t{}\t{}\t{}\t", std::str::from_utf8(self.header.tid2name(self.tid))?, self.pos + 1, char::from(ref_base), nbases + nins + ndel }
        if self.seq_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.seq_buf)?}
            self.seq_buf.clear();
        }

        if self.qual_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.qual_buf)?}
            self.qual_buf.clear();
        }

        print! {"\n"}

        Ok(())
    }

    pub fn set_pileup(&mut self) -> Result<bool, Error> {
        assert!(self.remove_buf.is_empty());
        let mut generated = false;

        let mut ndel @ mut nins @ mut nbases = 0;
        let ref_base = match &self.ref_seq {
            Some(seq) => seq[self.pos].to_ascii_uppercase(),
            None => b'N',
        };

        for (i, r) in self.rbuf.rbuf.iter_mut().enumerate() {
            if r.rec.reference_end() - 1 < self.pos as i64 {
                self.coverage -= 1;
                self.remove_buf.push_back(i);
                continue;
            }

            let mut ipos: i32 = -1;
            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32, &mut ipos);
            self.coverage += 1;

            match ret {
                Pileup::Op(Cigar::Match(_)) => {
                    get_base_pileup(
                        &r.cstate,
                        &r.rec,
                        ipos as u32,
                        self.pos,
                        &mut self.seq_buf,
                        &mut self.qual_buf,
                        ref_base,
                    );

                    nbases += 1;
                }

                Pileup::Op(Cigar::Ins(_)) => nins += 1,

                Pileup::Op(Cigar::Del(l)) => {
                    if ipos != -1 {
                        write_del(self.pos, &mut self.seq_buf, l as usize)?;
                        ndel += 1;
                    } else {
                        self.seq_buf.push(b'*');
                    }
                }

                Pileup::BeforePos() => {
                    panic!(
                        "{} {} {}",
                        r.rec.is_unmapped(),
                        self.pos,
                        r.rec.reference_end() - 1
                    );
                }

                Pileup::BaseEmpty() => (),
                _ => panic!(),
            }
        }

        if nbases + nins + ndel > 0 {
            self.write_pileup_str(ref_base, nbases, nins, ndel)?;
            generated = true;
        }

        while let Some(i) = self.remove_buf.pop_back() {
            self.rbuf.rbuf.swap_remove(i);
        }

        Ok(generated)
    }

    pub fn init_to_ref(&mut self) -> Result<IterResult, Error> {
        // todo: check if this works for bam files without refs in header
        //
        if self.tid == UNINIT_TID {
            self.tid = 0;
        } else {
            self.tid += 1;
        }

        if self.tid >= self.header.target_count() {
            Ok(IterResult::NoData)
        } else {
            self.max_pos = self.header.target_len(self.tid).context("No ref len")? as usize;
            self.pos = 0;
            self.next_pos = 0;
            self.reader.fetch((self.tid, 0, u32::MAX))?;
            Ok(IterResult::Generated)
        }
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        if self.pos >= self.max_pos {
            return Ok(IterResult::ReferenceEnd);
        }

        let mut gen = false;

        // if we are at the next position in the bam where reads are within window range,
        // resume read intake
        if self.pos == self.next_pos {
            let _r = self.fill_buffer();
        }

        // if we have reads in buffer, attempt to generate plp.
        if !self.rbuf.rbuf.is_empty() {
            gen = self.set_pileup()?;
        }

        // if no reads in buffer overlapped with pos, print empty plp if enabled
        if !gen && self.show_all {
            self.write_pileup_str(b'N', 0, 0, 0)?;
        }

        // if we need to print blank plps for each col,
        // advance query coord by 1
        // else, jump to the next coord with reads in range
        if self.show_all || !self.rbuf.rbuf.is_empty() {
            self.pos += 1;
        } else {
            self.pos = self.next_pos;
        }

        return Ok(IterResult::Generated);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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

        let mut ipos = 0;

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 1,
        };

        let mut ret = cigar_get_pos(&mut cstate, 4, &mut ipos);
        assert_eq!(ret, Pileup::Op(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, 5, &mut ipos);
        assert_eq!(ret, Pileup::Op(Cigar::Match(1)))
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

        let mut ipos = 0;

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 104,
        };

        let mut ret = cigar_get_pos(&mut cstate, 107, &mut ipos);
        assert_eq!(ret, Pileup::Op(Cigar::Match(4)));

        ret = cigar_get_pos(&mut cstate, 108, &mut ipos);
        assert_eq!(ret, Pileup::Op(Cigar::Ins(2)));

        ret = cigar_get_pos(&mut cstate, 109, &mut ipos);
        assert_eq!(ret, Pileup::Op(Cigar::Match(3)));
    }
}
