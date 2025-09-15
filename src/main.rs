use crate::params::parse_or_quit;
use crate::position_queue::{create_region_queue, PositionQueue};
use anyhow::Error;

mod alignment;
mod bamio;
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
mod utils;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    let mut pileup = pileup_iterator::PileupIterator::new(&params)?;

    if let Some(regstr) = params.inp.region {
        let reg_iter = create_region_queue(&regstr, &pileup.reader.header)?;
        pileup._auto_loop(&reg_iter)
    } else {
        let queue = PositionQueue::new(&pileup.reader.header)?;
        pileup._auto_loop(&queue)
    }
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
