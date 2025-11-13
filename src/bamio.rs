use crate::utils::has_index;

use anyhow::{Context, Error};
use rust_htslib::bam::{HeaderView, IndexedReader, Read, Reader, Record};
use std::path;

const READ_LENGTH_SAMPLE_SIZE: i8 = 10;

#[derive(Clone)]
pub enum BamDataSource {
    File(std::path::PathBuf),
    Stdin,
}

impl BamDataSource {
    pub fn from_string(s: &str) -> Result<Self, Error> {
        if s == "-" {
            Ok(Self::Stdin)
        } else if path::Path::exists(path::Path::new(s)) {
            Ok(Self::File(path::PathBuf::from(s)))
        } else {
            anyhow::bail!("Input path {} not found!", s)
        }
    }
}

impl std::fmt::Display for BamDataSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::File(f) => f.to_str().unwrap_or("FILE"),
            Self::Stdin => "STDIN",
        })
    }
}

pub struct BamReader {
    inner: Box<dyn BamRead>,
    pub header: HeaderView,
    pub cur_ref: String,
}

impl BamReader {
    pub fn new(src: &BamDataSource, threads: usize) -> Result<Self, Error> {
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

        Ok(Self { inner, header, cur_ref })
    }

    pub fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.inner.read_no_alloc(stored_read)
    }

    pub fn init_to_ref(&mut self, tid: u32, start: i64, end: i64) -> Result<(), Error> {
        self.cur_ref = std::str::from_utf8(self.header.tid2name(tid))?.to_string();
        self.inner.init_to_ref(tid, start, end)
    }

    pub fn sample_read_length(src: &BamDataSource) -> Result<usize, Error> {
        let mut temp_reader = Self::new(src, 1)?;

        let mut alloc = Record::new();

        let mut max_read_len: usize = 0;

        let mut reads_to_sample = READ_LENGTH_SAMPLE_SIZE;

        while reads_to_sample >= 0 {
            if let Some(r) = temp_reader.read_no_alloc(&mut alloc) {
                r?;
                max_read_len = std::cmp::max(max_read_len, alloc.seq_len());
                reads_to_sample -= 1;
            } else {
                break;
            }
        }

        if reads_to_sample == READ_LENGTH_SAMPLE_SIZE {
            anyhow::bail!("Failed to find any reads to sample for length! Is file {} empty?", src)
        }

        assert!(max_read_len > 0);
        Ok(max_read_len)
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
        self.read(stored_read).map(|e| e.context("Failed to retrieve read"))
    }
}

// Indexed Bam Reader
impl BamRead for IndexedReader {
    fn init_to_ref(&mut self, tid: u32, start: i64, end: i64) -> Result<(), Error> {
        self.fetch((tid, start, end)).context("Failed to fetch")
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
                anyhow::bail!("Attempted to create indexed reader from stdout!")
            }
        }
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read).map(|e| e.context("Failed to retrieve read"))
    }
}
