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

const UNINIT_POS: i64 = i64::MAX - 1;
const UNINIT_TID: i32 = i32::MAX - 1;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

pub struct PileupIterator {
    tid: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,
    tid_count: i32,
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

// cap qualitites at max of 126; this also helps avoid non-ascii output
pub fn get_qual(qual: u8) -> u8 {
    match qual.cmp(&92).is_gt() {
        true => 126,
        false => qual + 33,
    }
}

pub fn get_base_to_ref(
    cur_base: u8,
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
    pos: i64,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
) -> Result<(), Error> {
    let ipos = ipos as usize;
    let bam_pos = r.pos();

    if pos == bam_pos {
        seq_buf.push(FIRST_POS);
        seq_buf.push(get_qual(r.mapq()));
    }

    let cur_base = r.seq()[ipos];

    let base = get_base_to_ref(cur_base, pos as u64, refseq, r.is_reverse())?;

    let cur_qual = r.qual()[ipos];

    seq_buf.push(base);

    if pos == r.reference_end() - 1 {
        seq_buf.push(LAST_POS);
    }

    qual_buf.push(get_qual(cur_qual));

    Ok(())
}

pub fn write_del(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: i64,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
    del_len: i64,
) -> Result<(), Error> {
    write_match(cs, r, ipos, pos, seq_buf, qual_buf, refseq)?;
    write!(seq_buf, "-{}", del_len)?;
    for p in pos + 1..pos + del_len + 1 {
        let b = match refseq {
            Some(refseq) => get_base(refseq.get_base(p as u64)?, r.is_reverse()),
            None => b'N',
        };
        seq_buf.push(b);
    }
    Ok(())
}

pub fn write_ins(
    cs: &CigarState,
    r: &Record,
    ipos: u32,
    pos: i64,
    seq_buf: &mut Vec<u8>,
    qual_buf: &mut Vec<u8>,
    refseq: Option<&RefSeq>,
) -> Result<(), Error> {
    write_match(cs, r, ipos, pos, seq_buf, qual_buf, refseq)?;
    let mut k = cs.icig + 1;
    let ncig = cs.cig.len();
    while k < ncig {
        match cs.cig[k] {
            Cigar::Pad(l) => {
                seq_buf.extend(std::iter::repeat_n(b'*', l as usize));
            }

            Cigar::Ins(l) => {
                write!(seq_buf, "+{}", l)?;
                let (s, e) = ((ipos + 1) as usize, (ipos + 1 + l) as usize);
                for i in s..e {
                    let base = get_base(r.seq()[i], r.is_reverse());
                    seq_buf.push(base);
                }
            }

            _ => break,
        }

        k += 1;
    }

    Ok(())
}

pub fn cigar_get_pos(cs: &mut CigarState, pos: u32) -> CigarAtPos {
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

                cs.del = false;
                cs.qpos = pos as usize - cs.bam_pos as usize + cs.iseq as usize;
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
                cs.del = true;
                // cs.qpos = (cs.iseq + len) as usize;
                cs.qpos = cs.iseq as usize;
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
        let rbuf = read_buf::ReadBuffer::new(params.inp.depth, params.plp.disable_overlaps);
        let header = reader.header().clone();
        let show_all = params.plp.show_empty_coords;
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));
        let cur_rec = Record::new();
        let mut refseq = None;
        let min_baseq = params.plp.min_baseq;
        let max_tid = header.target_count() as i32;

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
            tid_count: max_tid,
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

    pub fn auto_loop(&mut self) -> Result<(), Error> {
        self.init_to_ref()?;

        loop {
            match self.next()? {
                IterResult::NoData => break,
                IterResult::Generated => continue,
                IterResult::ReferenceEnd => _ = self.init_to_ref()?,
            }
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
        print! {"{}\t{}\t{}\t{}\t", std::str::from_utf8(self.header.tid2name(self.tid as u32))?, self.pos + 1, char::from(ref_base), nbases + nins + ndel }
        if self.seq_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.seq_buf)?}
            self.seq_buf.clear();
        }

        if self.qual_buf.is_empty() {
            print! {"*"}
        } else {
            print! {"{}", std::str::from_utf8(&self.qual_buf)?}
            self.qual_buf.clear();
        }

        print! {"\n"}

        Ok(())
    }

    pub fn set_pileup(&mut self) -> Result<bool, Error> {
        assert!(self.rbuf.backup_buf.is_empty());
        let mut generated = false;
        let mut qual: u8;

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
                self.rbuf.depth -= 1;
                continue;
            }

            let ret = cigar_get_pos(&mut r.cstate, self.pos as u32);

            qual = if r.cstate.qpos >= r.rec.inner.core.l_qseq as usize {
                0
            } else {
                r.rec.qual()[r.cstate.qpos]
            };

            if qual < self.min_baseq {
                drop(r);
                self.rbuf.backup_buf.push(raw);
                continue;
            }

            // println! {"POS: {} QNAME: {}", self.pos, std::str::from_utf8(r.rec.qname())?}
            match ret {
                CigarAtPos::Op(Cigar::Match(_)) => {
                    write_match(
                        &r.cstate,
                        &r.rec,
                        r.cstate.qpos as u32,
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
                        r.cstate.qpos as u32,
                        self.pos,
                        &mut self.seq_buf,
                        &mut self.qual_buf,
                        self.refseq.as_ref(),
                    )?;
                }

                CigarAtPos::Op(Cigar::Del(l)) => {
                    if !r.cstate.del {
                        // write_del(self.pos, &mut self.seq_buf, l as usize)?;
                        write_del(
                            &r.cstate,
                            &r.rec,
                            r.cstate.qpos as u32,
                            self.pos,
                            &mut self.seq_buf,
                            &mut self.qual_buf,
                            self.refseq.as_ref(),
                            l as i64,
                        )?
                    } else {
                        self.seq_buf.push(b'*');
                        self.qual_buf.push(qual)
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

        let tlen = self.header.target_len(self.tid as u32).with_context(|| {
            format!(
                "Failed to get target length for {}",
                std::str::from_utf8(self.header.tid2name(self.tid as u32)).unwrap()
            )
        })?;

        if let Some(r) = self.refseq.as_mut() {
            // right now we just get the entire reference sequence.
            // Next step will be to load it in windows.
            let tidname = std::str::from_utf8(self.header.tid2name(self.tid as u32));
            r.load_seq(tidname?, 0, tlen)?
        }

        self.max_pos = tlen as i64;
        self.pos = 0;
        self.next_pos = 0;
        self.reader.fetch((self.tid, 0, u32::MAX))?;
        Ok(IterResult::Generated)
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        while let Some(read) = self.reader.read(&mut self.cur_rec) {
            read?;
            let r = &self.cur_rec;

            if r.is_unmapped() {
                continue;
            }

            if !self.read_filter.check_read(&r) {
                continue;
            }

            // TODO: resolve the int conversion mess
            if r.pos() < self.pos || r.tid() < self.tid {
                panic!("UNSORTED BAM")
            }

            let ret = self.rbuf.attempt_push(&r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::Unmapped => panic!(),
                read_buf::BufPushResult::DifferentReference => return Ok(IterResult::ReferenceEnd),
                read_buf::BufPushResult::MaxDepthMet => continue,
                read_buf::BufPushResult::Pushed => self.next_pos = r.pos(),
            }

            while self.pos < self.next_pos {
                self.set_pileup()?;
                self.pos += 1;
            }
        }

        while self.pos < self.max_pos {
            self.set_pileup()?;
            self.pos += 1;
        }

        match self.tid + 1 == self.tid_count as i32 {
            true => Ok(IterResult::NoData),
            false => Ok(IterResult::ReferenceEnd),
        }
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

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 1,
            qpos: 0,
            del: false,
        };

        let mut ret = cigar_get_pos(&mut cstate, 4);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, 5);
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

        let mut cstate = CigarState {
            cig: record.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: 104,
            qpos: 0,
            del: false,
        };

        let mut ret = cigar_get_pos(&mut cstate, 107);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(4)));

        ret = cigar_get_pos(&mut cstate, 108);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Ins(2)));

        ret = cigar_get_pos(&mut cstate, 109);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(3)));
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
        assert_eq!(ret, CigarAtPos::Op(Cigar::Del(4)));

        ret = cigar_get_pos(&mut cstate, 105);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 106);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 107);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 108);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Del(4)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 109);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(3)));
        assert_eq!(cstate.qpos, 1);

        ret = cigar_get_pos(&mut cstate, 110);
        assert_eq!(ret, CigarAtPos::Op(Cigar::Match(3)));
        assert_eq!(cstate.qpos, 2);
    }
}
