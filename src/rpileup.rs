use crate::read_buf;
use crate::read_buf::CigarState;
use anyhow::{Context, Error};
use rust_htslib::bam::record::{Cigar, Record};
use rust_htslib::bam::{HeaderView, Read, Reader};
use std::collections::HashSet;

pub struct PileUp {
    tid: u32,
    pos: usize,
    rbuf: read_buf::ReadBuffer,
    reader: Reader,
    header: HeaderView,
}

pub enum IterResult {
    ReferenceEnd,
    Generated,
    NoData,
}

#[derive(Debug, PartialEq, Eq)]
pub enum CigarResult {
    OutOfBounds(),
    Op(Cigar),
}

pub fn cigar_get_pos(cs: &mut CigarState, read: &Record, pos: u32, ipos: &mut i32) -> CigarResult {
    let cig = &cs.cig;
    let ncig = cig.len();
    // let cig = read.cigar();
    // let ncig = read.cigar_len();
    while cs.bam_pos < pos {
        if cs.icig >= ncig {
            return CigarResult::OutOfBounds();
        }

        let op = cig[cs.icig];
        match op {
            Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                let end_pos = cs.bam_pos + len + 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.iseq += len;
                    cs.icig += 1;
                }

                *ipos = pos as i32 - cs.bam_pos as i32 + cs.iseq as i32;
                if end_pos == pos && cs.icig + 1 < ncig {
                    let next_op = cig[cs.icig + 1];

                    match next_op {
                        Cigar::Ins(_) => return CigarResult::Op(next_op),
                        Cigar::Del(_) => return CigarResult::Op(next_op),
                        _ => (),
                    }
                }

                return CigarResult::Op(Cigar::Match(op.len()));
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

                *ipos -= 1;
                return CigarResult::Op(op);
            }

            Cigar::RefSkip(len) => {
                let end_pos = cs.bam_pos + len - 1;
                if end_pos < pos {
                    cs.bam_pos += len;
                    cs.icig += 1;
                    continue;
                }
            }
            _ => (),
        }
    }

    CigarResult::OutOfBounds()
}

impl PileUp {
    pub fn new(bam_fname: &str, tid: Option<usize>, pos: Option<usize>) -> Result<Self, Error> {
        let tid = tid.unwrap_or(0) as u32;
        let pos = pos.unwrap_or(0);
        let reader = Reader::from_path(bam_fname)?;
        let mut rbuf = read_buf::ReadBuffer::new();
        let header = reader.header().clone();
        rbuf.pos = pos;
        rbuf.tid = tid;

        Ok(Self {
            tid,
            pos,
            rbuf,
            reader,
            header,
        })
    }

    pub fn fill_buffer(&mut self) -> Result<read_buf::BufPushResult, Error> {
        let mut ret: read_buf::BufPushResult = read_buf::BufPushResult::DifferentReference;
        let mut pos: i64;

        if self.rbuf.rbuf.is_empty() {
            let first = self.reader.records().next().context("no read")??;
            self.rbuf.pos = first.pos() as usize;
            self.rbuf.tid = first.tid() as u32;
        }

        for rec in self.reader.records() {
            let r = rec?;

            ret = self.rbuf.push(r);

            match ret {
                read_buf::BufPushResult::AfterWindow => {
                    break;
                }
                read_buf::BufPushResult::DifferentReference => break,
                _ => (),
            }
        }
        Ok(ret)
    }

    pub fn set_pileup(&mut self) {
        let mut ndel @ mut nins @ mut nbases = 0;
        let mut to_remove: HashSet<usize> = HashSet::new();
        println! {"{}", self.rbuf.rbuf.len()};

        for (i, r) in self.rbuf.rbuf.iter_mut().enumerate() {
            let mut ipos: i32 = -1;
            let ret = cigar_get_pos(&mut r.cstate, &r.rec, self.pos as u32, &mut ipos);
            println! {"{:?}", ret}
            match ret {
                CigarResult::Op(Cigar::Match(_)) => {
                    let base = r.rec.seq().encoded_base(ipos as usize);
                    print! {" {base}"}
                    nbases += 1;
                }

                CigarResult::Op(Cigar::Ins(_)) => nins += 1,

                CigarResult::Op(Cigar::Del(_)) => {
                    if ipos != -1 {
                        ndel += 1;
                    }
                }

                CigarResult::OutOfBounds() => {
                    to_remove.insert(i);
                }

                _ => panic!(),
            }
        }

        print! {"\n"}
    }

    pub fn out_test(&self) {
        for r in self.rbuf.rbuf.iter() {
            println! {"{} {} {}", r.rec.pos(), r.rec.seq_len(), r.rec.cigar()};
        }
    }

    pub fn next(&mut self) -> Result<IterResult, Error> {
        self.pos += 1;
        if self.pos
            >= self
                .header
                .target_len(self.tid)
                .context("Unable to get ref len")? as usize
        {
            self.tid += 1;
            if self.header.target_count() <= self.tid {
                Ok(IterResult::NoData)
            } else {
                Ok(IterResult::ReferenceEnd)
            }
        } else {
            let _ = self.fill_buffer();
            // self.out_test();
            self.set_pileup();
            Ok(IterResult::Generated)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::CigarString;

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
            bam_pos: 0,
        };

        let mut ret = cigar_get_pos(&mut cstate, &record, 4, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Match(4)));
        ret = cigar_get_pos(&mut cstate, &record, 5, &mut ipos);
        assert_eq!(ret, CigarResult::Op(Cigar::Equal(1)))
    }
}
