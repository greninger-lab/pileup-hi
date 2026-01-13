use std::{fs::OpenOptions, io::BufWriter};

use crate::{
    alignment::PileupAlignment,
    bamio::OutputDataDest,
    engine::{MIN_BAM_READ_THREADS, MIN_COORDS_PER_THREAD},
};

pub type OutputWriter = BufWriter<Box<dyn std::io::Write>>;

use anyhow::{Context, Error};

pub fn read_ends_before_pos(a: &PileupAlignment, pos: i64) -> bool {
    a.rec.pos() + a.cstate.read_len_from_cigar - 1 < pos
}

pub fn temp_fname(prefix: &str, suffix: &str, ext: &str) -> String {
    format!("{prefix}_{suffix}.{ext}")
}

/// Get a writer to a particular destination. Lock specifies whether or not
/// we expect the writer to be the sole writer the source (pertinent if writing to stdout).
pub fn get_writer(handle: &OutputDataDest, writer_cap: usize, lock: bool, append: bool) -> Result<OutputWriter, Error> {
    let dest: Box<dyn std::io::Write> = match handle {
        OutputDataDest::File(p) => {
            let mut o = OpenOptions::new();
            Box::new(o.write(true).create(true).append(append).open(p)?)
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

pub struct FileThreadScheme {
    pub read_threads: usize,
    pub worker_threads: usize,
}

pub fn determine_thread_scheme(threads: usize, reflen: i64) -> FileThreadScheme {
    let reflen = reflen as usize;

    if reflen < MIN_COORDS_PER_THREAD {
        FileThreadScheme {
            read_threads: threads - 1,
            worker_threads: 1,
        }
    } else if reflen > MIN_COORDS_PER_THREAD * threads {
        FileThreadScheme {
            read_threads: MIN_BAM_READ_THREADS,
            worker_threads: threads,
        }
    } else {
        let mut nthreads = reflen / MIN_COORDS_PER_THREAD;
        let remainder = reflen % MIN_COORDS_PER_THREAD;

        if remainder >= MIN_COORDS_PER_THREAD {
            nthreads += 1;
        }

        FileThreadScheme {
            worker_threads: nthreads,
            read_threads: (threads - nthreads).max(MIN_BAM_READ_THREADS),
        }
    }
}
