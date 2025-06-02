use clap::Parser;

#[derive(Parser)]
pub struct Params {
    #[clap(flatten)]
    pub inp: InputParams,

    #[clap(flatten)]
    pub plp: PileupParams,

    #[clap(flatten)]
    pub outp: OutputParams,
}

#[derive(Parser)]
pub struct InputParams {
    #[arg(index = 1)]
    pub input: String,

    #[arg(short = 'f', long = "--fasta-ref")]
    pub refseq: Option<String>,

    #[arg(long = "tid")]
    pub tid: Option<u32>,

    #[arg(long = "pos")]
    pub pos: Option<usize>,
}

#[derive(Parser)]
pub struct PileupParams {
    #[arg(short = 'a')]
    pub show_empty_coords: bool,

    #[arg(long = "aa")]
    pub show_everything: bool,

    #[arg(long = "--rf")]
    pub incl_flags: Vec<String>,

    #[arg(long = "--ff", default_values_t = ["BAM_FSECONDARY".to_string(), "BAM_FQCFAIL".to_string(), "BAM_FDUP".to_string()])]
    pub excl_flags: Vec<String>,

    #[arg(short = 'Q', long = "--min-BQ", default_value_t = 13)]
    pub min_mapq: u8,
    // pub min_mapq: usize,
    // #[arg(short = 'd', default_value_t = 8000)]
    // pub max_depth: usize,
    // pub remove_overlaps: bool,
    // pub count_orphans: bool,
    // pub baq: bool,
}

#[derive(Parser)]
pub struct OutputParams {
    // pub output_ends: bool,
    // pub reverse_del: bool,
    // pub output_qname: bool,
}
