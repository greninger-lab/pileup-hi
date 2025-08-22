use crate::params::InputParams;
use crate::utils::has_index;

use anyhow::{Context, Error};
use rust_htslib::bam::{Format, Header, HeaderView, IndexedReader, Read, Reader, Record, Writer};

pub struct BamWriter {
    inner: Writer,
    _write_func: fn(&mut Self, &Record) -> Result<(), Error>,
}

// TODO: consider the overhead of using Option<Writer>. Wondering if it would be better to use
// function pointer to write nothing instead.
impl BamWriter {
    pub fn new_from_template(header: &HeaderView, output: &str) -> Result<Self, Error> {
        let header = Header::from_template(header);
        let inner = Writer::from_path(std::path::Path::new(output), &header, Format::Bam)?;
        let _write_func = Self::_write_read;

        Ok(Self { inner, _write_func })
    }

    pub fn write_record(&mut self, record: &Record) -> Result<(), Error> {
        // I have to do this wierd Ok() wrapping because Result<()> return type
        (self._write_func)(self, record)
    }

    pub fn void(header: &HeaderView) -> Result<Self, Error> {
        let header = Header::from_template(header);
        let inner = Writer::from_path(std::path::Path::new("/dev/null"), &header, Format::Sam)?;
        let _write_func = Self::_discard_read;

        Ok(Self { inner, _write_func })
    }

    fn _discard_read(&mut self, _rec: &Record) -> Result<(), Error> {
        Ok(())
    }
    fn _write_read(&mut self, rec: &Record) -> Result<(), Error> {
        Ok(self.inner.write(rec)?)
    }
}

pub struct BamReader {
    inner: Box<dyn BamRead>,
    pub header: HeaderView,
    pub max_tid: i32,
    pub cur_ref: String,
}

impl BamReader {
    pub fn new(params: &InputParams) -> Result<Self, Error> {
        match has_index(&params.input)? {
            true => {
                // println! {"Found index for {}.", &params.input}
                let inner = IndexedReader::new(&params.input, params.threads)?;
                let header = inner.header().clone();
                let max_tid = header.target_count() as i32;
                let cur_ref = "UNINIT".to_string();
                Ok(Self {
                    inner,
                    header,
                    max_tid,
                    cur_ref,
                })
            }

            false => {
                // println! {"No index found for {}. Using slower iteration...", &params.input}
                let inner = Reader::new(&params.input, params.threads)?;
                let header = inner.header().clone();
                let max_tid = header.target_count() as i32;
                let cur_ref = "UNINIT".to_string();
                Ok(Self {
                    inner,
                    header,
                    max_tid,
                    cur_ref,
                })
            }
        }
    }

    pub fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.inner.read_no_alloc(stored_read)
    }

    pub fn init_to_ref(&mut self, tid: u32) -> Result<(), Error> {
        self.cur_ref = std::str::from_utf8(self.header.tid2name(tid))?.to_string();
        self.inner.init_to_ref(tid)
    }
}

/// An interface used to allow reading both indexed and un-indexed bams with the same struct.
pub trait BamRead {
    fn init_to_ref(&mut self, tid: u32) -> Result<(), Error>;
    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized;
    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>>;
}

// Standard BAM Reader NO INDEX
impl BamRead for Reader {
    fn init_to_ref(&mut self, _tid: u32) -> Result<(), Error> {
        Ok(())
    }

    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret = Self::from_path(input_file)?;
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read)
            .map(|e| e.context("Failed to retrieve read"))
    }
}

// Indexed Bam Reader
impl BamRead for IndexedReader {
    fn init_to_ref(&mut self, tid: u32) -> Result<(), Error> {
        self.fetch((tid, 0, u32::MAX)).context("Failed to fetch")
    }

    fn new(input_file: &str, threads: usize) -> Result<Box<Self>, Error>
    where
        Self: Sized,
    {
        let mut ret = Self::from_path(input_file)?;
        ret.set_threads(threads)?;
        Ok(Box::new(ret))
    }

    fn read_no_alloc(&mut self, stored_read: &mut Record) -> Option<Result<(), Error>> {
        self.read(stored_read)
            .map(|e| e.context("Failed to retrieve read"))
    }
}
