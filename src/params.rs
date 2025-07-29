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

pub fn validate_params(p: &Params) {
    if p.plp.indel_realign && p.inp.refseq.is_none() {
        eprintln!{"You must provide a reference fasta when enabling indel realignment."};
        std::process::exit(1);
    }
}

pub fn parse_or_quit() -> Params {
    match Params::try_parse() {
        Ok(p) => {validate_params(&p); p},
        Err(e) => {
            e.print().unwrap();
            std::process::exit(1)
        }
    }
}

#[derive(Parser)]
pub struct InputParams {
    #[arg(index = 1)]
    pub input: String,

    #[arg(short = 't', long = "--threads", default_value_t = num_cpus::get())]
    pub threads: usize,

    #[arg(short = 'f', long = "fasta-ref")]
    pub refseq: Option<String>,

    #[arg(long = "tid")]
    pub tid: Option<i32>,

    #[arg(long = "pos")]
    pub pos: Option<i64>,

    #[arg(short = 'd', long = "depth", default_value_t = 8000)]
    pub depth: usize,
}

#[derive(Parser)]
pub struct PileupParams {
    #[arg(short = 'a')]
    pub show_empty_coords: bool,

    #[arg(short = 'i')]
    pub indel_realign: bool,

    #[arg(long = "aa")]
    pub show_everything: bool,

    #[arg(short = 'x')]
    pub disable_overlaps: bool,

    #[arg(short = 'A')]
    pub count_orphans: bool,

    #[arg(long = "rf")]
    pub incl_flags: Vec<String>,

    #[arg(long = "ff", default_values_t = ["BAM_FSECONDARY".to_string(), "BAM_FQCFAIL".to_string(), "BAM_FDUP".to_string(), "BAM_FUNMAP".to_string()])]
    pub excl_flags: Vec<String>,

    #[arg(short = 'q', long = "min-MQ", default_value_t = 0)]
    pub min_mapq: u8,

    #[arg(short = 'Q', long = "min-BQ", default_value_t = 13)]
    pub min_baseq: u8,

    #[arg(short = 'x', long = "disable_overlap_removal")]
    pub disable_overlap: bool,
    // pub baq: bool,
}

#[derive(Parser)]
pub struct OutputParams {
    // pub output_ends: bool,
    // pub reverse_del: bool,
    // pub output_qname: bool,
}
