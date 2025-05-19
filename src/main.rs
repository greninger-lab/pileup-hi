use crate::read_buf::{BufPushResult, ReadBuffer};
use anyhow::Error;
use clap::Parser;
use rust_htslib::bam::{Header, Read, Reader, Record, Writer};

mod read_buf;
mod rpileup;

#[derive(Parser)]
pub struct Args {
    pub input: String,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let pos = 782;
    let tid = 0;

    let mut pileup = rpileup::PileUp::new(&args.input, Some(tid), Some(pos))?;
    pileup.next()?;

    Ok(())
}
