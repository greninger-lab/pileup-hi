use crate::params::parse_or_quit;
use anyhow::Error;

mod bamio;
mod overlap;
mod params;
mod pileup;
mod pileup_iterator;
mod pileup_writer;
mod read_buf;
mod read_filter;
mod read_walker;
mod realigner;
mod refseq;
mod utils;
mod left_align_indel;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    let mut pileup = pileup_iterator::PileupIterator::new(params)?;
    pileup.auto_loop()
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
