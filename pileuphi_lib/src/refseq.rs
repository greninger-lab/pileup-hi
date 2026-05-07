use crate::errors::{Error, ErrorKind};
use bio::io::{
    fasta,
    fasta::{FastaRead, Record},
};
use log::warn;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::sync::{Arc, Weak};

struct FastaIndexedReader {
    inner: fasta::IndexedReader<File>,
}

struct FastaReader {
    inner: fasta::Reader<BufReader<File>>,
}

pub type RefSeqHandle = Arc<Option<Vec<u8>>>;

pub trait ReadsFasta {
    fn read_to_bytes(&mut self, refname: &str) -> Result<Option<Vec<u8>>, Error>;
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
                warn!("Unable to find ref {refname} in fasta. Proceeding without reference...");
                return Ok(None); // read through all refs without finding one matching the given id
            }
        }
    }
}

//////////////////////////////////////////
/// Holds reference information requested by any number of processing threads, giving threads
/// read-only access to a reference on demand. Responsible for freeing unused references and loading
/// new ones.
//////////////////////////////////////////
pub struct RefSeq {
    data: RefCell<HashMap<String, Weak<Option<Vec<u8>>>>>,
    file: String,
}

// internally, reference sequences are held by atomic reference-counted pointers. If we don't have
// any threads using the reference anymore, we free the bytes storing the sequence. This is only
// really efficient at preventing constant re-loading when we make our threads process interval
// jobs in reference-sorted order (which we do unless we are given a BED file or some other custom
// list of intervals).

impl RefSeq {
    pub fn get_reader(file: &str) -> Result<Box<dyn ReadsFasta>, Error> {
        // TODO: use regular fasta reader to avoid using rust_bio
        let idx_name = format! {"{file}.fai"};
        let faidx = Path::new(&idx_name);

        let reader: Box<dyn ReadsFasta> = if !faidx.exists() {
            Box::new(FastaReader {
                inner: fasta::Reader::from_file(Path::new(&file)).map_err(|e| {
                    Error::from_generic(e.into(), ErrorKind::RefSeqError("Failed to create FASTA reader"))
                })?,
            })
        } else {
            Box::new(FastaIndexedReader {
                inner: fasta::IndexedReader::from_file(&Path::new(&file)).map_err(|e| {
                    Error::from_generic(
                        e.into(),
                        ErrorKind::RefSeqError("Failed to create indexed FASTA reader"),
                    )
                })?,
            })
        };
        Ok(reader)
    }

    pub fn load_seq(file_name: &str, ref_name: &str) -> Result<Option<Vec<u8>>, Error> {
        let mut reader = RefSeq::get_reader(file_name)?;
        reader.read_to_bytes(ref_name)
    }

    pub fn new(file: String) -> Self {
        Self {
            data: RefCell::new(HashMap::new()),
            file,
        }
    }

    pub fn yield_handle(&self, ref_name: &str) -> Result<RefSeqHandle, Error> {
        let mut lock = self.data.borrow_mut();

        if let Some(slot) = lock.get_mut(ref_name) {
            if let Some(copy) = slot.upgrade() {
                Ok(copy)
            } else {
                let data = Arc::new(RefSeq::load_seq(&self.file, ref_name)?);
                *slot = Arc::downgrade(&data);
                Ok(data)
            }
        } else {
            let data = Arc::new(RefSeq::load_seq(&self.file, ref_name)?);
            lock.insert(ref_name.to_string(), Arc::downgrade(&data));
            Ok(data)
        }
    }
}
