use anyhow::Error;
use clap::Parser;

mod read_buf;
mod rpileup;

#[derive(Parser)]
pub struct Args {
    pub input: String,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();
    let pos = 0;
    let tid = 0;

    let mut pileup = rpileup::PileUp::new(&args.input, Some(tid), Some(pos))?;
    let mut ret: rpileup::IterResult;
    loop {
        ret = pileup.next()?;
        match ret {
            rpileup::IterResult::NoData => break,
            _ => (),
        }
    }

    Ok(())
}
