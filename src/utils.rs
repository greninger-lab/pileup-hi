use std::{fs::OpenOptions, io::BufWriter};

use crate::{alignment::PileupAlignment, bamio::OutputDataDest};

pub type OutputWriter = BufWriter<Box<dyn std::io::Write>>;

use anyhow::{Context, Error};

pub fn read_ends_before_pos(a: &PileupAlignment, pos: i64) -> bool {
    a.rec.pos() + a.cstate.read_len_from_cigar - 1 < pos
}

pub fn temp_fname(prefix: &str, suffix: &str, ext: &str) -> String {
    format!("{prefix}_{suffix}.{ext}")
}

/// Get a writer to a particular destination. Lock specifies whether or not
/// we expect the writer to be the sole writer the source
pub fn get_writer_multi(
    handle: &OutputDataDest,
    writer_cap: usize,
    lock: bool,
    append: bool,
) -> Result<OutputWriter, Error> {
    let dest: Box<dyn std::io::Write> = match handle {
        OutputDataDest::File(p) => {
            let mut o = OpenOptions::new();
            let file = o.write(true).create(true).append(append).open(p)?;

            if lock {
                file.lock()?;
            }
            Box::new(file)
        }

        OutputDataDest::Stdout => {
            if lock {
                Box::new(std::io::stdout().lock())
            } else {
                Box::new(std::io::stdout())
            }
        }
    };

    Ok(BufWriter::with_capacity(writer_cap, dest))
}

pub fn has_index(bam_file: &str) -> Result<bool, Error> {
    let potential_index = format! {"{bam_file}.bai"};

    std::fs::exists(&potential_index)
        .with_context(|| format!("Unable to check directory for file {}", &potential_index))
}
