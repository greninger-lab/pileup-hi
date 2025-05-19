use crate::read_buf;
use anyhow::{Context, Error};
use rust_htslib::bam::{HeaderView, Read, Reader};

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

        for rec in self.reader.records() {
            let r = rec?;

            pos = r.pos();
            ret = self.rbuf.push(r);

            match ret {
                read_buf::BufPushResult::AfterWindow => {
                    println! {"{} {} {}", self.pos, self.rbuf.len, pos}
                    break;
                }
                read_buf::BufPushResult::DifferentReference => break,
                _ => (),
            }
        }
        Ok(ret)
    }

    pub fn out_test(&self) {
        for r in self.rbuf.rbuf.iter() {
            println! {"{} {} {}", r.pos(), r.seq_len(), r.cigar()};
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
            self.out_test();
            Ok(IterResult::Generated)
        }
    }
}
