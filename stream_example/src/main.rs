#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

mod args;

use pileuphi_lib::{outputs::PileupString, PileupCoordinate, PileupEngine};

use crate::args::parse_or_quit;

fn main() {
    let params = parse_or_quit();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let mut engine = PileupEngine::init_stream(params.plp, PileupString::new()).expect("failed to init pileup engine");

    let iters = engine.get_iter(params.inp).expect("Failed to get iterator");

    for mut iter in iters {
        while let Some(b) = iter.advance() {
            match b.expect("Error during pileup generation") {
                PileupCoordinate::NoCoverage => (),
                PileupCoordinate::Coverage(b) => b.write(&mut std::io::stdout()).unwrap(),
            }
        }
    }
}
