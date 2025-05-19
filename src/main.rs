use anyhow::Error;
use clap::Parser;
use rust_htslib::bam::{Header, Read, Reader, Record, Writer};

#[derive(Parser)]
pub struct Args {
    pub input: String,
    // pub output: Option<String>,
}

pub struct ReadBuffer {
    rbuf: Vec<Record>,
    len: usize,
    pos: usize,
    tid: u32,
}

impl ReadBuffer {
    pub fn push(&mut self, r: Record) -> u8 {
        if r.tid() as u32 != self.tid {
            0
        } else {
            if r.seq_len() > self.len {
                self.len = r.seq_len();
            }
            if r.pos() as usize + self.len < self.pos {
                1
            } else {
                if r.pos() as usize > self.pos + self.len {
                    1
                } else {
                    self.rbuf.push(r);
                    0
                }
            }
        }
    }

    pub fn new() -> Self {
        let rbuf: Vec<Record> = Vec::with_capacity(500);
        let len = 0;
        let pos = 0;
        let tid = 0;

        Self {
            rbuf,
            len,
            pos,
            tid,
        }
    }
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let mut reader = Reader::from_path(args.input)?;
    let pos = 782;
    let tid = 0;
    let mut read_buf = ReadBuffer::new();
    read_buf.tid = tid;
    read_buf.pos = pos;

    for record in reader.records() {
        read_buf.push(record?);
    }

    let mut writer = Writer::from_stdout(
        &Header::from_template(reader.header()),
        rust_htslib::bam::Format::Bam,
    )?;

    for r in read_buf.rbuf {
        // writer.write(&r)?
        println! {"POS: {}, LEN: {}", r.pos(), r.seq_len()};
    }

    Ok(())
}
