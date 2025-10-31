#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::params::parse_or_quit;
use crate::position_queue::PositionQueue;
use crate::threading::PileupMultiThreader;
use anyhow::Error;
use pileup_iterator::PileupIterator;

mod alignment;
mod bamio;
mod output;
mod overlap;
mod params;
mod pileup_iterator;
mod pileup_writer;
mod position_queue;
mod read_buf;
mod read_filter;
mod read_walker;
mod realigner;
mod refseq;
mod threading;
mod utils;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();

    let queue = PositionQueue::new_from_bam(&params.inp.input)?;

    if params.inp.threads > 1 {
        let mut driver = PileupMultiThreader::new(queue, params)?;
        driver.run()?;
    } else {
        eprintln!("Running in single-thread mode...");
        let mut iterator = PileupIterator::new(&params, None)?;
        iterator._auto_loop(&queue)?;
    }

    Ok(())
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
