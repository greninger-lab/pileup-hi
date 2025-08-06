use clap::Parser;
use indel_left_align::{bam_from_file, left_align_indels, load_seq, make_bam_out};
use rust_htslib::bam::{record::CigarString, Read, Record};

#[derive(Parser)]
pub struct Args {
    #[arg(index = 1)]
    pub input: String,

    #[arg(index = 2)]
    pub output: String,

    #[arg(index = 3)]
    pub reference: String,
}

pub fn main() {
    let args = Args::parse();

    let mut reader = bam_from_file(&args.input);
    let mut writer = make_bam_out(&args.output, &reader);
    let mut cig: CigarString;
    let mut rec = Record::new();

    let refseq = load_seq(&args.reference, None);

    while let Some(r) = reader.read(&mut rec) {
        r.unwrap();
        cig = left_align_indels(&rec, refseq.as_slice());
        rec.set_cigar(Some(&cig));
        writer.write(&rec).unwrap();
    }
}
