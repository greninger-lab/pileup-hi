use crate::alignment::PileupAlignment;
use anyhow::{Context, Error};

pub fn read_ends_before_pos(a: &PileupAlignment, pos: i64) -> bool {
    a.rec.pos() + a.cstate.read_len_from_cigar - 1 < pos
}

pub fn temp_fname(prefix: &str, suffix: &str, ext: &str) -> String {
    format!("{prefix}_{suffix}.{ext}")
}

pub fn has_index(bam_file: &str) -> Result<bool, Error> {
    let potential_index = format! {"{bam_file}.bai"};

    std::fs::exists(&potential_index)
        .with_context(|| format!("Unable to check directory for file {}", &potential_index))
}
