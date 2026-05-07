#[cfg(feature = "cli")]
use {crate::engine::MIN_COORDS_PER_THREAD, crate::read_filter::BamFlag, clap::Parser};

pub const STDOUT_ARG_STR: &str = "STDOUT";

#[cfg_attr(feature = "cli", derive(Parser, Clone))]
pub struct InputParams {
    /// File to read (path) or stdout ("-")
    pub file: String,

    /// only process a particular bam region (e.g. chr1:0-8000)
    #[cfg_attr(feature = "cli", arg(short = 'r', long = "region"))]
    pub region: Option<String>,
}

#[cfg_attr(feature = "cli", derive(Parser))]
#[derive(Clone)]
pub struct PileupParams {
    #[cfg_attr(feature = "cli", arg(short = 'a'))]
    pub show_empty_coords: bool,

    /// Output positions for regions with no depth
    #[cfg_attr(feature = "cli", arg(long = "aa"))]
    pub show_empty_regions: bool,

    #[cfg_attr(feature = "cli", arg(short = 'c', long = "thread-coords", default_value_t = MIN_COORDS_PER_THREAD))]
    pub coords_per_thread: i64,

    /// Reference fasta to use for comparison, must be indexed
    #[cfg_attr(feature = "cli", arg(short = 'f', long = "fasta-ref"))]
    pub refseq: Option<String>,

    /// The maximum number of reads to sample per position. Set to 0 to uncap
    #[cfg_attr(feature = "cli", arg(short = 'd', long = "depth", default_value_t = 8000))]
    pub depth: usize,

    /// Disable R1/R2 mate overlap correction
    #[cfg_attr(feature = "cli", arg(short = 'x'))]
    pub disable_overlaps: bool,

    /// Count reads with unmapped mates
    #[cfg_attr(feature = "cli", arg(short = 'A'))]
    pub count_orphans: bool,

    #[cfg_attr(feature = "cli", arg(long = "rf"))]
    pub incl_flags: Vec<BamFlag>,

    /// Don't consider any reads with these flags
    #[cfg_attr(feature = "cli", arg(long = "ff", default_values_t = [BamFlag::Secondary, BamFlag::QCFail, BamFlag::Duplicate, BamFlag::Unmapped]))]
    pub excl_flags: Vec<BamFlag>,

    /// Minimum mapping quality for a read's bases to be counted
    #[cfg_attr(feature = "cli", arg(short = 'q', long = "min-MQ", default_value_t = 0))]
    pub min_mapq: u8,

    /// Minimum phred score for a base to be counted
    #[cfg_attr(feature = "cli", arg(short = 'Q', long = "min-BQ", default_value_t = 13))]
    pub min_baseq: u8,

    /// Disable calcluation of base alignment quality (BAQ)
    #[cfg_attr(feature = "cli", arg(short = 'B', long = "no-BAQ", default_value_t = false))]
    pub no_baq: bool,

    /// Calculate BAQ even when BAQ already exists
    #[cfg_attr(
        feature = "cli",
        arg(short = 'E', long = "redo-BAQ", default_value_t = false, conflicts_with("no_baq"))
    )]
    pub redo_baq: bool,
}
