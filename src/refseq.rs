use anyhow::{Context, Error};
use bio::io::{
    fasta,
    fasta::{FastaRead, Record},
};
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

// This is necessary because rust-bio's fasta reader reads to different types depending on whether
// the fasta is indexed. I need to change this in the future.
pub enum RefSeqStore {
    Record(Record),
    Seq(Vec<u8>),
}
pub trait FastaReader {
    fn read_to_bytes(&mut self, refname: &str, seqbuf: &mut RefSeqStore) -> Result<i64, Error>;
}

impl FastaReader for fasta::Reader<BufReader<File>> {
    fn read_to_bytes(&mut self, refname: &str, seqbuf: &mut RefSeqStore) -> Result<i64, Error> {
        match seqbuf {
            RefSeqStore::Record(r) => {
                while let Ok(()) = self.read(r) {
                    if r.id() == refname {
                        return Ok(r.seq().len() as i64);
                    }

                    if r.seq().is_empty() {
                        anyhow::bail!("Unable to find sequence {} in file", refname);
                    }
                }

                anyhow::bail!("Unable to find sequence {} in file", refname);
            }

            RefSeqStore::Seq(_) => {
                anyhow::bail!("Cannot read from non-indexed reader into plain byte slice!")
            }
        }
    }
}

impl FastaReader for fasta::IndexedReader<File> {
    fn read_to_bytes(&mut self, refname: &str, seqbuf: &mut RefSeqStore) -> Result<i64, Error> {
        match seqbuf {
            RefSeqStore::Record(_) => {
                anyhow::bail!("Cannot read from indexed reader into Record struct")
            }
            RefSeqStore::Seq(seq) => {
                self.fetch_all(refname)?;
                self.read(seq)
                    .with_context(|| format!("Fail to read reference {}", refname))?;
                Ok(seq.len() as i64)
            }
        }
    }
}

pub struct RefSeq {
    wend: i64,
    store: RefSeqStore,
    reader: Box<dyn FastaReader>,
}

impl RefSeq {
    // TODO: use regular fasta reader to avoid using rust_bio
    pub fn from_file(file: &str) -> Result<Self, Error> {
        // check if idx exists
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        if !faidx.exists() {
            let reader = Box::new(fasta::Reader::from_file(Path::new(&file))?);
            let store = RefSeqStore::Record(Record::new());
            Ok(Self { wend: 0, store, reader })
        } else {
            let reader = Box::new(fasta::IndexedReader::from_file(&Path::new(&file))?);
            let store = RefSeqStore::Seq(Vec::new());
            Ok(Self { wend: 0, store, reader })
        }
    }

    pub fn load_seq(&mut self, t_name: &str) -> Result<(), Error> {
        let reflen = self.reader.read_to_bytes(t_name, &mut self.store)?;
        self.wend = reflen;
        Ok(())
    }

    #[allow(dead_code)]
    pub fn get_base(&self, pos: i64) -> Result<u8, Error> {
        if pos > self.wend {
            anyhow::bail!("Position {pos} exceeds current loaded window up to {}", self.wend)
        } else {
            match &self.store {
                RefSeqStore::Record(r) => Ok(r.seq()[pos as usize]),

                RefSeqStore::Seq(seq) => Ok(seq[pos as usize]),
            }
        }
    }

    pub fn len(&self) -> i64 {
        self.wend
    }

    pub fn yield_seq(&self) -> &[u8] {
        match &self.store {
            RefSeqStore::Record(r) => r.seq(),
            RefSeqStore::Seq(seq) => seq,
        }
    }
}
