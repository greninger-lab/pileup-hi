use anyhow::{Context, Error};
use bio::io::fasta;
use std::fs::File;
use std::path::Path;

pub struct RefSeq {
    wstart: u64,
    wend: u64,
    seq: Vec<u8>,
    // file: String,
    reader: fasta::IndexedReader<File>,
}

impl RefSeq {
    pub fn from_file(file: String) -> Result<Self, Error> {
        // check if idx exists
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        if !faidx.exists() {
            anyhow::bail!(
                "Unable to find index file {:?} for ref file {file}",
                idx_name
            )
        }

        let reader = fasta::IndexedReader::from_file(&Path::new(&file))?;
        let wstart = 0;
        let wend = 0;
        let seq = Vec::new();

        Ok(Self {
            wstart,
            wend,
            seq,
            reader,
        })
    }

    pub fn load_seq(&mut self, t_name: &str, start: u64, stop: u64) -> Result<(), Error> {
        self.reader.fetch(t_name, start, stop).with_context(|| {
            format!(
                "Internal error: unable to fetch interval {} - {} for ref {}",
                start, stop, t_name,
            )
        })?;

        self.wstart = start;
        self.wend = stop;

        self.reader.read(&mut self.seq).with_context(|| {
            format!(
                "ref fasta corruption at interval {} - {} for ref {}",
                start, stop, t_name
            )
        })?;

        Ok(())
    }

    pub fn get_base(&self, pos: u64) -> Result<u8, Error> {
        if self.seq.is_empty() {
            return Ok(b'x');
        }

        if pos > self.wend {
            anyhow::bail!(
                "Position {pos} exceeds current loaded window up to {}",
                self.wend
            )
        } else {
            let offset = pos - self.wstart;
            let r = self
                .seq
                .get(offset as usize)
                .context("Unable to get ref base")?;
            Ok(*r)
        }
    }
}
