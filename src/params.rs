use clap::{Parser, Subcommand};

#[derive(Subcommand, Clone)]
/// Tuple of command + recommended threads
pub enum Commands {
    /// Generate a samtools mpileup string
    Plp(Params),
    /// Generate a per-coordinate count of bases and indels
    Histo(Params),
}

#[derive(Parser, Clone)]
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

    /// Override for # of processing threads, changing can result in slowdown
    #[arg(short = 't', long = "threads")]
    pub threads: Option<usize>,
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
    #[arg(short = 'a')]
    pub show_empty_coords: bool,

    /// Reference fasta to use for comparison, must be indexed
    #[arg(short = 'f', long = "fasta-ref")]
    pub refseq: Option<String>,

    /// The maximum number of reads to sample per position. Set to 0 to uncap
    #[arg(short = 'd', long = "depth", default_value_t = 8000)]
    pub depth: usize,

    /// Output positions for regions with no depth
    #[arg(long = "aa")]
    pub show_everything: bool,

    /// Disable R1/R2 mate overlap correction
    #[arg(short = 'x')]
    pub disable_overlaps: bool,

    /// Count reads with unmapped mates
    #[arg(short = 'A')]
    pub count_orphans: bool,

    #[arg(long = "rf")]
    pub incl_flags: Vec<String>,

    /// Don't any reads with these flags
    #[arg(long = "ff", default_values_t = ["BAM_FSECONDARY".to_string(), "BAM_FQCFAIL".to_string(), "BAM_FDUP".to_string(), "BAM_FUNMAP".to_string()])]
    pub excl_flags: Vec<String>,

    /// Minimum mapping quality for a read's bases to be counted
    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    pub min_mapq: u8,

    /// Minimum phred score for a base to be counted
    #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
    pub min_baseq: u8,
}
