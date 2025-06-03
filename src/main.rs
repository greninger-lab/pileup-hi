use crate::params::parse_or_quit;
use anyhow::Error;

mod overlap;
mod params;
mod pileup;
mod read_buf;
mod read_filter;
mod refseq;
mod rpileup;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();

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

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
