use crate::params::Params;
use crate::pileup::CigarState;
use crate::read_buf;
use crate::read_filter::ReadFilter;
use crate::refseq::RefSeq;
use anyhow::{Context, Error};
use num_cpus;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{ext::BamRecordExtensions, HeaderView, IndexedReader, Read, Record};
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
    refseq: Option<RefSeq>,
    seq_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    read_filter: ReadFilter,
    cur_rec: Record,
    min_baseq: u8,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarAtPos {
    BeforePos(),
    Op(Cigar),
    BaseEmpty(),
}

pub fn get_base(mut cur_base: u8, is_reverse: bool) -> u8 {
    match is_reverse {
        false => cur_base.make_ascii_uppercase(),
        true => cur_base.make_ascii_lowercase(),
    }

    cur_base
}

pub fn get_base_to_ref(
    mut cur_base: u8,
    ref_coord: u64,
    refseq: Option<&RefSeq>,
    is_reverse: bool,
) -> Result<u8, Error> {
    if let Some(refseq) = refseq {
        let ref_base = refseq.get_base(ref_coord)?;
        if ref_base == cur_base {
            if is_reverse {
                Ok(R_MATCH)
            } else {
                Ok(F_MATCH)
            }
        } else {
            Ok(get_base(cur_base, is_reverse))
        }
    } else {
        Ok(get_base(cur_base, is_reverse))
    }
}

pub fn write_match(
    _cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
) -> Result<(), Error> {
    // assert_ne!(ipos, -1);
    let ipos = ipos as usize;
    // let bam_pos = cs.bam_pos as usize;
    let bam_pos = r.pos() as usize;

    if pos == bam_pos {
        seq_buf.push(FIRST_POS);
        seq_buf.push(r.mapq() + 33);
    }

    let cur_base = r.seq()[ipos];

    let base = get_base_to_ref(cur_base, pos as u64, refseq, r.is_reverse())?;

    let cur_qual = r.qual()[ipos] + 33;

    seq_buf.push(base);

    if pos == r.reference_end() as usize - 1 {
        seq_buf.push(LAST_POS);
    }

    qual_buf.push(cur_qual);

    Ok(())
}

pub fn write_del(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
    del_len: usize,
) -> Result<(), Error> {
    write_match(cs, r, ipos, pos, seq_buf, qual_buf, refseq)?;
    write!(seq_buf, "-{}", del_len)?;
    for _ in pos..pos + del_len {
        seq_buf.push(b'N')
    }
    Ok(())
}

pub fn write_ins(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: usize,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
) -> Result<(), Error> {
    write_match(cs, r, ipos, pos, seq_buf, qual_buf, refseq)?;
    let mut k = cs.icig + 1;
    let ncig = cs.cig.len();
    let ipos = ipos + 1;
    while k < ncig {
        match cs.cig[k] {
            Cigar::Pad(l) => {
                seq_buf.extend(std::iter::repeat_n(b'*', l as usize));
            }

            Cigar::Ins(l) => {
                write!(seq_buf, "+{}", l)?;
                let (s, e) = (ipos as usize, (ipos + l) as usize);
                for i in s..e {
                    let base = get_base(r.seq()[i], r.is_reverse());
                    seq_buf.push(base);
                    qual_buf.push(r.qual()[i] + 33);
                }
            }

            _ => break,
        }

        k += 1;
    }

    Ok(())
}

pub fn cigar_get_pos(cs: &mut CigarState, pos: u32, ipos: &mut i32) -> CigarAtPos {
    let cig = &cs.cig;
    let ncig = cig.len();
    let mut op: Cigar;
    while cs.bam_pos <= pos {
        if cs.icig >= ncig {
            // this should never happen, since we check cigars beforehand to at least end
            // at the queried coordinate, if not pass over it.
            return CigarAtPos::BeforePos();
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
                        Cigar::Ins(_) => return CigarAtPos::Op(next_op),
                        Cigar::Del(_) => return CigarAtPos::Op(next_op),
                        _ => (),
                    }
                }
                return CigarAtPos::Op(Cigar::Match(len));
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
                return CigarAtPos::Op(op);
            }

            Cigar::RefSkip(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }

                return CigarAtPos::BaseEmpty();
            }
            _ => (),
        }
    }

    CigarAtPos::BaseEmpty()
}

impl PileupIterator {
    pub fn new(params: Params) -> Result<Self, Error> {
        let tid = params.inp.tid.unwrap_or(UNINIT_TID);
        let pos @ next_pos @ max_pos = params.inp.pos.unwrap_or(UNINIT_POS);
        let mut reader = IndexedReader::from_path(params.inp.input)?;
        reader.set_threads(num_cpus::get())?;
        let rbuf = read_buf::ReadBuffer::new();
        let header = reader.header().clone();
        let show_all = params.plp.show_empty_coords;
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));
        let cur_rec = Record::new();
        let mut refseq = None;
        let min_baseq = params.plp.min_baseq;

        let read_filter = ReadFilter::new(
            params.plp.min_mapq,
            params.plp.count_orphans,
            params.plp.excl_flags.iter().map(|s| s.as_str()).collect(),
            params.plp.incl_flags.iter().map(|s| s.as_str()).collect(),
        )?;

        if let Some(ref_file) = params.inp.refseq {
            refseq = Some(RefSeq::from_file(ref_file)?);
        }

        Ok(Self {
            tid,
            pos,
            next_pos,
            max_pos,
            rbuf,
            reader,
            header,
            min_baseq,
            read_filter,
            show_all,
            refseq,
            seq_buf,
            qual_buf,
            cur_rec,
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
        // let mut scanned = 0;

        let mut prev_pos = i64::MIN;

        if self.cur_rec.tid() == -1 {
            if let Some(rec) = self.reader.read(&mut self.cur_rec) {
                rec?;
            } else {
                // if we have no reads at all to set next pos, assume
                // we've hit the end of the reference, and set next pos to MAX
                self.next_pos = usize::MAX;
                return Ok(());
            }
        }

        loop {
            if self.cur_rec.tid() == -1 {
                break;
            }

            let r = &self.cur_rec;

            if r.is_unmapped() {
                match self.reader.read(&mut self.cur_rec) {
                    None => {
                        self.cur_rec.set_tid(-1);
                        break;
                    }
                    Some(_) => continue,
                };
            }

            if !self.read_filter.check_read(&r) {
                match self.reader.read(&mut self.cur_rec) {
                    None => {
                        self.cur_rec.set_tid(-1);
                        break;
                    }

                    Some(_) => continue,
                };
            }

            if r.pos() < prev_pos {
                panic!("UNSORTED BAM! {} {}", r.pos(), prev_pos)
            }

            prev_pos = r.pos();

            ret = self.rbuf.attempt_push(&r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::Unmapped => panic!(),

                read_buf::BufPushResult::AfterWindow(next_pos) => {
                    self.next_pos = self.pos + next_pos;
                    break;
                }

                read_buf::BufPushResult::DifferentReference => {
                    break;
                }
                read_buf::BufPushResult::Pushed => {
                    self.cur_rec.set_tid(-1);
                    match self.reader.read(&mut self.cur_rec) {
                        Some(Ok(_)) => continue,
                        None => {
                            self.cur_rec.set_tid(-1);
                            break;
                        }
                        Some(Err(_)) => panic!(),
                    }
                }
            };
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
        assert!(self.rbuf.backup_buf.is_empty());
        let mut generated = false;

        let mut ndel @ mut nins @ mut nbases = 0;
        let ref_base = match &self.refseq {
            Some(seq) => seq.get_base(self.pos as u64)?,
            None => b'N',
        };

        for raw in self.rbuf.rbuf.drain(..) {
            let mut r = raw.borrow_mut();

            if r.rec.reference_end() - 1 < self.pos as i64 {
                drop(r);
                drop(raw);
                continue;
            }

            let mut ipos: i32 = 0;
            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32, &mut ipos);

            let qual_idx = if ipos == -1 {
                if r.cstate.iseq == 0 {
                    drop(r);
                    self.rbuf.backup_buf.push(raw);
                    continue;
                } else {
                    r.cstate.iseq as usize
                }
            } else {
                ipos as usize
            };

            if r.rec.qual()[qual_idx] < self.min_baseq {
                drop(r);
                self.rbuf.backup_buf.push(raw);
                continue;
            }

            match ret {
                CigarAtPos::Op(Cigar::Match(_)) => {
                    write_match(
                        &r.cstate,
                        &r.rec,
                        qual_idx as u32,
                        self.pos,
                        &mut self.seq_buf,
                        &mut self.qual_buf,
                        self.refseq.as_ref(),
                    )?;

                    nbases += 1;
                }

                CigarAtPos::Op(Cigar::Ins(_)) => {
                    nins += 1;
                    write_ins(
                        &r.cstate,
                        &r.rec,
                        // ipos as u32,
                        qual_idx as u32,
                        self.pos,
                        &mut self.seq_buf,
                        &mut self.qual_buf,
                        self.refseq.as_ref(),
                    )?;
                }

                CigarAtPos::Op(Cigar::Del(l)) => {
                    if ipos != -1 {
                        // write_del(self.pos, &mut self.seq_buf, l as usize)?;
                        write_del(
                            &r.cstate,
                            &r.rec,
                            // ipos as u32,
                            qual_idx as u32,
                            self.pos,
                            &mut self.seq_buf,
                            &mut self.qual_buf,
                            self.refseq.as_ref(),
                            l as usize,
                        )?
                    } else {
                        self.seq_buf.push(b'*');
                    }
                    ndel += 1;
                }

                CigarAtPos::BeforePos() => {
                    panic!(
                        "{} {} {}",
                        r.rec.is_unmapped(),
                        self.pos,
                        r.rec.reference_end() - 1
                    );
                }

                CigarAtPos::BaseEmpty() => (),
                _ => panic!(),
            }

            drop(r);
            self.rbuf.backup_buf.push(raw);
        }

        if nbases + nins + ndel > 0 {
            self.write_pileup_str(ref_base, nbases, nins, ndel)?;
            generated = true;
        }

        self.rbuf.reset();

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
            if let Some(r) = self.refseq.as_mut() {
                // right now we just get the entire reference sequence.
                // Next step will be to load it in windows.
                let tidname = std::str::from_utf8(self.header.tid2name(self.tid));
                r.load_seq(
                    tidname?,
                    0,
                    self.header
                        .target_len(self.tid)
                        .context("Failed to get target length")?,
                )?
            }

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
        // if self.pos == self.next_pos {
        //     self.fill_buffer()?;
        // }

        self.fill_buffer()?;

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
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, 5, &mut ipos);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(1)))
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
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(4)));

        ret = cigar_get_pos(&mut cstate, 108, &mut ipos);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Ins(2)));

        ret = cigar_get_pos(&mut cstate, 109, &mut ipos);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(3)));
    }
}
