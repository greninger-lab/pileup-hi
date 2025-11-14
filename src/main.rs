#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::params::parse_or_quit;
use crate::threading::PileupEngine;
use anyhow::Error;

mod alignment;
mod bamio;
mod basedepth_string;
mod cigar_resolve;
mod output;
mod overlap;
mod params;
mod pileup_iterator;
mod pileup_string;
mod position_queue;
mod read_buf;
mod read_filter;
mod read_walker;
mod refseq;
mod threading;
mod utils;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();

    let mut engine = PileupEngine::initialize(params.inp, params.plp)?;
    engine.run()?;

    Ok(())
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
