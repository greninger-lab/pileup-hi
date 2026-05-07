use crate::{params::STDOUT_ARG_STR, utils::has_index};

use crate::errors::{Error, ErrorKind};
use rust_htslib::bam::{HeaderView, IndexedReader, Read, Reader, Record};
use std::path;

const BAM_READ_LEN_SAMPLE_SIZE: usize = 10;

#[derive(Clone)]
pub enum BamDataSource {
    File(std::path::PathBuf),
    Stdin,
}

#[derive(Clone)]
pub enum OutputDataDest {
    File(String),
    Stdout,
}

impl BamDataSource {
    pub fn has_index(&self) -> Result<bool, Error> {
        match self {
            Self::File(f) => has_index(f.to_str().unwrap_or("")),
            Self::Stdin => Ok(false),
        }
    }

    // get everything before extension
    pub fn fname(&self) -> Result<String, Error> {
        match self {
            Self::File(f) => {
                let full = f.to_str().unwrap_or("");
                if let Some((fname, _)) = full.rsplit_once(".") {
                    Ok(fname.to_string())
                } else {
                    Ok(full.to_string())
                }
            }

            Self::Stdin => Ok(STDOUT_ARG_STR.to_string()),
        }
    }
}

impl BamDataSource {
    pub fn from_string(s: &str) -> Result<Self, Error> {
        if s == "-" {
            Ok(Self::Stdin)
        } else if path::Path::exists(path::Path::new(s)) {
            Ok(Self::File(path::PathBuf::from(s)))
        } else {
            Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                std::io::ErrorKind::NotFound,
                "BAM not found",
            ))))
        }
    }
}

impl std::fmt::Display for BamDataSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::File(f) => f.to_str().unwrap_or("FILE NAME INVALID (character encoding?)"),
            Self::Stdin => "STDIN",
        })
    }
}

impl OutputDataDest {
    pub fn from_string(s: &str) -> Self {
        if s == STDOUT_ARG_STR {
            Self::Stdout
        } else {
            Self::File(s.to_string())
        }
    }
}

impl std::fmt::Display for OutputDataDest {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::File(f) => f,
            Self::Stdout => "STDOUT",
        })
    }
}

pub struct BamReader {
    inner: Box<dyn BamRead>,
    has_index: bool,
    pub header: HeaderView,
    pub cur_ref: String,
    pub eof: bool,
    pub src: BamDataSource,
}

impl BamReader {
    pub fn new(src: &BamDataSource, threads: usize) -> Result<Self, Error> {
        let _has_index = src.has_index()?;

        let inner: Box<dyn BamRead> = match &src {
            BamDataSource::File(file) => match has_index(file.to_str().unwrap())? {
                true => {
                    let mut inner = IndexedReader::new(src, threads)?;
                    inner.fetch(".")?;
                    inner
                }

                false => Reader::new(src, threads)?,
            },

            BamDataSource::Stdin => Reader::new(src, threads)?,
        };

        let header = inner.get_header().clone();
        let cur_ref = "UNINIT".to_string();

        Ok(Self {
            inner,
            header,
            cur_ref,
            eof: false,
            has_index: _has_index,
            src: src.clone(),
        })
    }

    pub fn sample_read_len(src: &BamDataSource) -> Result<usize, Error> {
        let mut reader = Self::new(src, 1)?; // 1 thread
        let mut cached = Record::new();

        let mut nsampled = 0;
        let mut totallen = 0;

        for _ in 0..BAM_READ_LEN_SAMPLE_SIZE {
            if let Some(r) = reader.read_no_alloc(&mut cached) {
                r?;
                nsampled += 1;
                totallen += cached.seq_len();
            } else {
                break;
            }
        }

        if nsampled == 0 {
            Ok(0)
        } else {
            Ok(totallen / nsampled)
        }
    }

    pub fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        // if we call read_no_alloc() on an unindexed reader after it already returned None (eof),
        // it will infinitely hang, at least with the version of rust-htslib I'm using.
        if self.eof {
            return None;
        }

        self.inner.read_no_alloc(stored_read)
    }

    pub fn init_to_ref(&mut self, tid: u32, start: i64, end: i64) -> Result<(), Error> {
        self.cur_ref = std::str::from_utf8(self.header.tid2name(tid))?.to_string();
        self.inner.init_to_ref(tid, start, end)?;

        if self.eof && self.has_index {
            self.eof = false;
        }

        Ok(())
    }
}

/// An interface used to allow reading both indexed and un-indexed bams with the same struct.
pub trait BamRead {
    fn init_to_ref(&mut self, tid: u32, start: i64, end: i64) -> Result<(), Error>;
    fn get_header(&self) -> &HeaderView;
    fn new(src: &BamDataSource, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized;
    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>>;
}

// Standard BAM Reader NO INDEX
impl BamRead for Reader {
    fn init_to_ref(&mut self, _tid: u32, _start: i64, _end: i64) -> Result<(), Error> {
        Ok(())
    }

    fn get_header(&self) -> &HeaderView {
        self.header()
    }

    fn new(src: &BamDataSource, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret;
        match src {
            BamDataSource::File(f) => {
                ret = Self::from_path(f)?;
            }
            BamDataSource::Stdin => {
                ret = Self::from_stdin()?;
            }
        }
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read).map(|res| res.map_err(Error::from))
    }
}

// Indexed Bam Reader
impl BamRead for IndexedReader {
    fn init_to_ref(&mut self, tid: u32, start: i64, end: i64) -> Result<(), Error> {
        self.fetch((tid, start, end)).map_err(Error::from)
    }

    fn get_header(&self) -> &HeaderView {
        self.header()
    }

    fn new(src: &BamDataSource, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret;
        match src {
            BamDataSource::File(f) => {
                ret = Self::from_path(f)?;
            }

            BamDataSource::Stdin => {
                return Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                    std::io::ErrorKind::InvalidFilename,
                    "Attemtped to create an indexed reader for STDIN",
                ))));
            }
        }
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read).map(|res| res.map_err(Error::from))
    }
}
