use anyhow::{Context, Error};
use bio::io::fasta;
use std::fs::File;
use std::path::Path;

pub struct RefSeq {
    wstart: u64,
    wend: u64,
    seq: Vec<u8>,
    reader: Option<fasta::IndexedReader<File>>,
    empty: bool,
}

impl RefSeq {
    pub fn is_empty(&self) -> bool {
        self.empty
    }

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

        let reader = Some(fasta::IndexedReader::from_file(&Path::new(&file))?);
        let wstart = 0;
        let wend = 0;
        let seq = Vec::new();

        Ok(Self {
            wstart,
            wend,
            seq,
            reader,
            empty: false,
        })
    }

    pub fn empty() -> Self {
        Self {
            wstart: 0,
            wend: 0,
            seq: vec![],
            reader: None,
            empty: true,
        }
    }

    pub fn load_seq(&mut self, t_name: &str, start: u64, stop: u64) -> Result<(), Error> {
        if let Some(reader) = &mut self.reader {
            reader.fetch(t_name, start, stop).with_context(|| {
                format!(
                    "Internal error: unable to fetch interval {} - {} for ref {}",
                    start, stop, t_name,
                )
            })?;

            self.wstart = start;
            self.wend = stop;

            reader.read(&mut self.seq).with_context(|| {
                format!(
                    "ref fasta corruption at interval {} - {} for ref {}",
                    start, stop, t_name
                )
            })?;

            Ok(())
        } else {
            anyhow::bail!("Attempted to fetch sequence with no reader loaded")
        }
    }

    pub fn get_base(&self, pos: u64) -> Result<u8, Error> {
        if self.empty {
            return Ok(b'N');
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

    pub fn get_interval(&mut self, start: u64, end: u64) -> Result<&[u8], Error> {
        if self.empty {
            self.seq = std::iter::repeat_n(b'N', (end - start + 1) as usize).collect();
            return Ok(self.seq.as_slice());
        }

        if start < self.wstart || end > self.wend {
            anyhow::bail!(
                "Invalid query window for loaded ref seq window {}-{}: {}-{}",
                self.wstart,
                self.wend,
                start,
                end
            )
        } else {
            Ok(&self.seq[start as usize..=end as usize])
        }
    }
}
