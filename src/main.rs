#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::{
    basedepth_string::BaseDepthString,
    output::setup_exit_handler,
    params::{parse_or_quit, Commands},
    pileup_string::PileupString,
    threading::PileupEngine,
};

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

const PLP_RECOMMENDED_THREADS: usize = 2;
const HISTO_RECOMMENDED_THREADS: usize = 4;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    setup_exit_handler();

    match params.command {
        Commands::Plp(params) => {
            let threads = params.threads.unwrap_or(PLP_RECOMMENDED_THREADS);
            let engine =
                PileupEngine::initialize(params.inp, params.plp, threads, PileupString::new())?;
            engine.run()?
        }

        Commands::Histo(params) => {
            let threads = params.threads.unwrap_or(HISTO_RECOMMENDED_THREADS);
            let engine =
                PileupEngine::initialize(params.inp, params.plp, threads, BaseDepthString::new())?;
            engine.run()?;
        }
    };

    Ok(())
}

fn main() {
    if let Err(e) = _main() {
        eprintln!("Error: {}", e);
        std::process::exit(1);
    }
}
