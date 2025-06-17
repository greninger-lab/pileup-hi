use crate::params::parse_or_quit;
use anyhow::Error;

mod overlap;
mod params;
mod pileup;
mod read_buf;
mod read_filter;
mod read_walker;
mod refseq;
mod rpileup;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    let mut pileup = rpileup::PileupIterator::new(params)?;
    pileup.auto_loop()
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
