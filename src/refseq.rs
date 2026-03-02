use anyhow::Error;
use bio::io::{
    fasta,
    fasta::{FastaRead, Record},
};
use log::warn;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::Arc;

struct FastaIndexedReader {
    inner: fasta::IndexedReader<File>,
}

struct FastaReader {
    inner: fasta::Reader<BufReader<File>>,
}

struct NoReader {}

impl NoReader {
    fn new() -> Self {
        Self {}
    }
}

pub trait ReadsFasta {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error>;
}

impl ReadsFasta for NoReader {
    fn read_to_bytes(&mut self, _refname: &str) -> Result<Option<Vec<u8>>, Error> {
        Ok(None)
    }
}

impl ReadsFasta for FastaIndexedReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error> {
        // ref not found
        if self.inner.fetch_all(refname).is_err() {
            warn!("Unable to find ref {refname} in fasta. Proceeding without reference...");
            return Ok(None);
        };

        let mut output = Vec::new();
        self.inner.read(&mut output)?;

        Ok(Some(output))
    }
}

impl ReadsFasta for FastaReader {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error> {
        let mut record: Record = Default::default();

        loop {
            self.inner.read(&mut record)?;
            if record.id() == refname {
                return Ok(Some(record.seq().to_vec())); // found it
            } else if record.seq().is_empty() {
                return Ok(None); // read through all refs without finding one matching the given id
            }
        }
    }
}

pub struct RefSeq {
    reader: Box<dyn ReadsFasta>,
    data: Option<Arc<Vec<u8>>>,
}

impl RefSeq {
    // TODO: use regular fasta reader to avoid using rust_bio
    pub fn from_file(file: &str) -> Result<Self, Error> {
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        if !faidx.exists() {
            let reader = FastaReader {
                inner: fasta::Reader::from_file(Path::new(&file))?,
            };

            Ok(Self {
                reader: Box::new(reader),
                data: None,
            })
        } else {
            let reader = FastaIndexedReader {
                inner: fasta::IndexedReader::from_file(&Path::new(&file))?,
            };

            Ok(Self {
                reader: Box::new(reader),
                data: None,
            })
        }
    }

    pub fn blank() -> Self {
        Self {
            reader: Box::new(NoReader::new()),
            data: None,
        }
    }

    pub fn load_seq(&mut self, t_name: &str) -> Result<(), Error> {
        self.data = if let Some(bytes) = self.reader.read_to_bytes(t_name)? {
            Some(Arc::new(bytes))
        } else {
            None
        };
        Ok(())
    }

    pub fn yield_handle(&self) -> Option<Arc<Vec<u8>>> {
        if let Some(ref refseq) = self.data {
            Some(Arc::clone(refseq))
        } else {
            None
        }
    }
}
