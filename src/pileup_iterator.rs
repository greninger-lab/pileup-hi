use crate::bamio::BamReader;
use crate::params::Params;
use crate::pileup::CigarState;
use crate::pileup_writer::PileupWriter;
use crate::read_buf;
use crate::read_filter::ReadFilter;
use crate::refseq::RefSeq;

use anyhow::{Context, Error};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};

const UNINIT_POS: i64 = i64::MAX - 1;
const UNINIT_TID: i32 = i32::MAX - 1;

/// The iterator responsible for coordinate-wise traversal of a BAM reference/region.
/// Maintains a buffer of pileups, iterator state (current reference and position), and
/// auxiliary structs needed for variant detection and output.
pub struct PileupIterator {
    tid: i32,
    pos: i64,
    next_pos: i64,
    max_pos: i64,
    tid_count: i32,
    show_all: bool,
    rbuf: read_buf::ReadBuffer,
    pileup_writer: PileupWriter,
    reader: BamReader,
    refseq: Option<RefSeq>,
    read_filter: ReadFilter,
    cur_rec: Record,
    min_baseq: u8,
}

/// A counter class to track the number of reads differing from the reference at a given pileup.
/// Holds the qualities of all bases, and the number of reads with ref-matching vs differing bases.
pub struct PileupContext {
    quals: Vec<u8>,
    nref: u32,
    nalt: u32,
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

/// Get the cigar operation in a read at a given index. Intended to mimic cigar_resolver2 from
/// htslib.
///
/// If the queried index is at the end of a match operation, the function will check if the next
/// operation is a deletion or insertion, and return the corresponding operation if so.
///
/// For example:
///
/// if return == [CigarAtPos(Cigar::Del(l))], then current position is [Cigar::Match] but the very next
/// one is [Cigar::Del].
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
        let reader = BamReader::new(&params.inp)?;
        let pileup_writer = PileupWriter::new();
        let rbuf = read_buf::ReadBuffer::new(params.inp.depth, params.plp.disable_overlaps);
        let show_all = params.plp.show_empty_coords;
        let cur_rec = Record::new();
        let mut refseq = None;
        let min_baseq = params.plp.min_baseq;
        let max_tid = reader.max_tid;

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
            pileup_writer,
            reader,
            min_baseq,
            read_filter,
            show_all,
            refseq,
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

    /// Generate a pileup from all bases passing the minimum quality filter at a given reference
    /// position.
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

            match ret {
                CigarAtPos::Op(Cigar::Match(_)) => {
                    self.pileup_writer.write_match(
                        &r.rec,
                        r.cstate.qpos as u32,
                        self.pos,
                        &self.refseq,
                    )?;

                    nbases += 1;
                }

                CigarAtPos::Op(Cigar::Ins(_)) => {
                    self.pileup_writer.write_match(
                        &r.rec,
                        r.cstate.qpos as u32,
                        self.pos,
                        &self.refseq,
                    )?;

                    self.pileup_writer
                        .write_insertion(&r.cstate, &r.rec, r.cstate.qpos as u32)?;
                    nins += 1;
                }

                CigarAtPos::Op(Cigar::Del(l)) => {
                    if !r.cstate.del {
                        self.pileup_writer.write_match(
                            &r.rec,
                            r.cstate.qpos as u32,
                            self.pos,
                            &self.refseq,
                        )?;

                        self.pileup_writer.write_deletion_start(
                            &r.rec,
                            self.pos + 1,
                            &self.refseq,
                            l as i64,
                        )?
                    } else {
                        self.pileup_writer.write_deletion(qual);
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
            self.pileup_writer.write_pileup_str(
                ref_base,
                self.pos,
                nbases,
                nins,
                ndel,
                &self.reader.cur_ref,
            )?;
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

        self.reader.init_to_ref(self.tid as u32)?;

        let tlen = self
            .reader
            .header
            .target_len(self.tid as u32)
            .with_context(|| format!("Failed to get target length for {}", self.reader.cur_ref))?;

        if let Some(r) = self.refseq.as_mut() {
            // right now we just get the entire reference sequence.
            // Next step will be to load it in windows.
            r.load_seq(&self.reader.cur_ref, 0, tlen)?
        }

        self.max_pos = tlen as i64;
        self.pos = 0;
        self.next_pos = 0;

        Ok(IterResult::Generated)
    }

    /// Iterate over reads until there is either no more input data left or a read mapped to a
    /// different reference is encountered.
    pub fn next(&mut self) -> Result<IterResult, Error> {
        while let Some(read) = self.reader.read_no_alloc(&mut self.cur_rec) {
            read?;
            let r = &self.cur_rec;

            if r.is_unmapped() {
                continue;
            }

            if !self.read_filter.check_read(&r) {
                continue;
            }

            if r.pos() < self.pos || r.tid() < self.tid {
                panic!("UNSORTED BAM")
            }

            let ret = self.rbuf.attempt_push(&r, self.pos, self.tid);

            match ret {
                read_buf::BufPushResult::Unmapped => panic!(),
                read_buf::BufPushResult::DifferentReference => return Ok(IterResult::ReferenceEnd),

                // if we hit depth limit, we exhaust all reads at this position to avoid dealing
                // with them at self.pos + 1.
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
