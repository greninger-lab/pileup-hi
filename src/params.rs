use clap::{crate_authors, crate_description, crate_version, Parser, Subcommand};
pub const STDOUT_ARG_STR: &str = "STDOUT";
use crate::engine::MIN_COORDS_PER_THREAD;

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

#[derive(Parser, Clone)]
pub struct Params {
    #[clap(flatten)]
    pub inp: InputParams,

    #[clap(flatten)]
    pub plp: PileupParams,
}

#[derive(Parser, Clone)]
pub struct InputParams {
    /// File to read (path) or stdout ("-")
    pub file: String,

    /// only process a particular bam region (e.g. chr1:0-8000)
    #[arg(short = 'r', long = "region")]
    pub region: Option<String>,
}

#[derive(Parser, Clone)]
pub struct PileupParams {
    /// Where to write all output to
    #[arg(short = 'o', long = "output", default_value_t = STDOUT_ARG_STR.to_string())]
    pub output: String,

    #[arg(short = 'a')]
    pub show_empty_coords: bool,

    /// Output positions for regions with no depth
    #[arg(long = "aa")]
    pub show_empty_regions: bool,

    /// Number of threads per reference
    #[arg(short = 't', long = "threads", default_value_t = 3)]
    pub threads: usize,

    #[arg(short = 'c', long = "thread-coords", default_value_t = MIN_COORDS_PER_THREAD)]
    pub coords_per_thread: i64,

    /// Reference fasta to use for comparison, must be indexed
    #[arg(short = 'f', long = "fasta-ref")]
    pub refseq: Option<String>,

    /// The maximum number of reads to sample per position. Set to 0 to uncap
    #[arg(short = 'd', long = "depth", default_value_t = 8000)]
    pub depth: usize,

    /// Disable R1/R2 mate overlap correction
    #[arg(short = 'x')]
    pub disable_overlaps: bool,

    /// Count reads with unmapped mates
    #[arg(short = 'A')]
    pub count_orphans: bool,

    #[arg(long = "rf")]
    pub incl_flags: Vec<String>,

    /// Don't consider any reads with these flags
    #[arg(long = "ff", default_values_t = ["BAM_FSECONDARY".to_string(), "BAM_FQCFAIL".to_string(), "BAM_FDUP".to_string(), "BAM_FUNMAP".to_string()])]
    pub excl_flags: Vec<String>,

    /// Minimum mapping quality for a read's bases to be counted
    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum phred score for a base to be counted
    #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
    pub min_baseq: u8,

    /// Disable calcluation of base alignment quality (BAQ)
    #[arg(short = 'B', long = "no-BAQ", default_value_t = false)]
    pub no_baq: bool,

    /// Calculate BAQ even when BAQ already exists
    #[arg(short = 'E', long = "redo-BAQ", default_value_t = false, conflicts_with("no_baq"))]
    pub redo_baq: bool,
}
