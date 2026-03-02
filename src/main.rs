#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

use crate::{
    basedepth_string::BaseDepthString,
    engine::PileupEngine,
    output::setup_exit_handler,
    params::{parse_or_quit, Commands},
    pileup_string::PileupString,
};

use anyhow::Error;
use log::error;

#[cfg(debug_assertions)]
use log::warn;

mod alignment;
mod bamio;
mod baq;
mod basedepth_string;
mod cigar_resolve;
mod engine;
mod output;
mod overlap;
mod params;
mod pileup_iterator;
mod pileup_string;
mod position_queue;
mod read_buf;
mod read_filter;
// mod realign_report;
mod refseq;
mod utils;

fn _main() -> Result<(), Error> {
    let params = parse_or_quit();
    setup_exit_handler();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    #[cfg(debug_assertions)]
    warn!("Running the debug build of pileup-hi! This is slow.");

    match params.command {
        Commands::Plp(params) => {
            let engine = PileupEngine::initialize(params.inp, params.plp, PileupString::new())?;
            engine.run()?
        }

        Commands::Histo(params) => {
            let engine = PileupEngine::initialize(params.inp, params.plp, BaseDepthString::new())?;
            engine.run()?;
        }
    };

    Ok(())
}

fn main() {
    if let Err(e) = _main() {
        error!("Error: {}", e);
        std::process::exit(1);
    }
}
