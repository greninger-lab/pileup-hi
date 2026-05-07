use clap::{crate_authors, crate_description, crate_version, Parser, Subcommand};

use pileuphi_lib::param::{InputParams, PileupParams, STDOUT_ARG_STR};

#[derive(Subcommand, Clone)]
pub enum Commands {
    /// Generate a samtools mpileup string
    Plp(Params),
    /// Generate a per-coordinate count of bases and indels
    Histo(Params),
}

#[derive(Parser, Clone)]
#[command(
    name = "pileup-hi",
    version = crate_version!(),
    author = crate_authors!(),
    about = crate_description!(),
    help_template = "===== {name} {version} ===== \n{about}\n{author}\n\n{usage-heading} {usage}\n\n{all-args}"
)]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Parser, Clone)]
pub struct Params {
    #[clap(flatten)]
    pub inp: InputParams,

    #[clap(flatten)]
    pub plp: PileupParams,

    /// Output to write to
    #[arg(short = 'o', long = "output", default_value_t = STDOUT_ARG_STR.to_string())]
    pub output: String,

    /// Number of threads per reference
    #[arg(short = 't', long = "threads", default_value_t = 3)]
    pub threads: usize,
}

pub fn parse_or_quit() -> Args {
    match Args::try_parse() {
        Ok(p) => {
            // no argument checking at the moment, leaving here for the future.
            p
        }
        Err(e) => {
            e.print().unwrap();
            std::process::exit(1)
        }
    }
}
