use crate::params::Params;
use anyhow::Error;
use clap::Parser;

mod params;
mod pileup;
mod read_buf;
mod read_filter;
mod refseq;
mod rpileup;

fn main() -> Result<(), Error> {
    let params = Params::parse();

    let mut pileup = rpileup::PileupIterator::new(params)?;
    let mut ret: rpileup::IterResult;

    loop {
        ret = pileup.init_to_ref()?;

        match ret {
            rpileup::IterResult::NoData => break,
            _ => loop {
                match pileup.next()? {
                    rpileup::IterResult::ReferenceEnd => break,
                    rpileup::IterResult::NoData => panic!(),
                    _ => (),
                }
            },
        }
    }

    Ok(())
}
