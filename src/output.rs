#![allow(dead_code)]
use crate::alignment::PileupAlignment;
use anyhow::Error;
use log::warn;
use std::fs::File;
use std::io::{stdout, BufReader, BufWriter, Write};
use std::sync::Mutex;

const PILEUP_OUTPUT_BUF_PURGE_THRES: usize = 1_000_000;

pub static TEMP_FILES: Mutex<Vec<String>> = Mutex::new(Vec::new());

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Get the reference of the pileup
    fn tid(&self) -> i32;
    /// Get the coordinate of the pileup
    fn pos(&self) -> i64;
    /// Update internal data with pileup alignment
    fn intake(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error>;
    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, ref_seq: Option<&[u8]>);
    fn write<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), Error>;
    fn depth(&self) -> u32;
    fn clear(&mut self);
    fn new() -> Self;
}

#[derive(Clone)]
pub struct PileupOutputChunk {
    pub data: Vec<u8>,
    pub id: usize,
}

pub struct TempOutputHandle {
    fname: String,
    writer: BufWriter<Box<dyn Write>>,
}

impl TempOutputHandle {
    pub fn new(fname: &str, writer_cap: usize) -> Result<Self, Error> {
        let dest: Box<dyn Write> = if fname == "STDOUT" {
            Box::new(stdout().lock())
        } else {
            Box::new(File::create(fname)?)
        };

        Ok(Self {
            fname: fname.to_string(),
            writer: BufWriter::with_capacity(writer_cap, dest),
        })
    }

    pub fn write(&mut self, data: &[u8]) {
        let _ = self.writer.write_all(data);
        self.writer.flush().unwrap();
    }

    pub fn delete(&mut self) -> Result<(), Error> {
        std::fs::remove_file(std::path::Path::new(&self.fname)).map_err(Error::msg)
    }
}

pub fn merge_temp_outputs<W: std::io::Write>(outputs: &[String], mut dest: W) -> Result<(), Error> {
    for fname in outputs {
        if fname == "STDOUT" {
            continue;
        }
        let fhandle = File::open(fname)?;
        let mut reader = BufReader::with_capacity(2 * 1024 * 1024, fhandle);
        std::io::copy(&mut reader, &mut dest)?;
        std::fs::remove_file(fname)?;
    }
    dest.flush()?;
    Ok(())
}

pub fn setup_exit_handler() {
    ctrlc::set_handler(|| {
        warn!("Received termination signal. Cleaning up intermediate files...");
        if let Ok(files) = TEMP_FILES.lock() {
            for fname in files.iter() {
                if fname == "STDOUT" {
                    continue;
                };
                let _ = std::fs::remove_file(fname);
            }
        }

        std::process::exit(130);
    })
    .expect("Failed to set exit handler")
}

/// A chunked dynamic array used for batching data writes and reducing allocations, intended for
/// multithreading where a worker thread also owns its writer. Chunks span a range of coordinates,
/// each of which should be assigned its output or None if the coordinate failed to meet an output
/// criterion (e.g. depth).
pub struct PileupOutputArray<T: OrderedPileupOutput> {
    data: Vec<Vec<Option<T>>>,
    capacity: usize,
    pub cur_entry: usize,
    cur_chunk: usize,
    write_batch_size: usize,
    output: TempOutputHandle,
    pub id: usize,
    outbuf: Vec<u8>,
}

impl<T: OrderedPileupOutput> PileupOutputArray<T> {
    pub fn alloc_chunk(&mut self) {
        let n_chunks = self.capacity / self.write_batch_size;
        let remainder = self.capacity % self.write_batch_size;

        self.data = Vec::with_capacity(n_chunks);

        for _ in 0..n_chunks - 1 {
            self.data.push(vec![Some(T::new()); self.write_batch_size]);
        }

        let final_size = remainder + self.write_batch_size;

        self.data.push(vec![Some(T::new()); final_size]);

        self.cur_entry = 0;
        self.cur_chunk = 0;
    }

    pub fn new(
        capacity: usize,
        write_batch_size: usize,
        id: usize,
        output: TempOutputHandle,
    ) -> Result<Self, Error> {
        let outbuf = Vec::with_capacity(write_batch_size * size_of::<T>());
        let mut s = Self {
            data: Vec::new(),
            capacity,
            cur_entry: 0,
            cur_chunk: 0,
            output,
            write_batch_size,
            outbuf,
            id,
        };

        s.alloc_chunk();
        Ok(s)
    }

    pub fn get_current_mut(&mut self) -> &mut T {
        self.data[self.cur_chunk][self.cur_entry].as_mut().unwrap()
    }

    pub fn tombstone(&mut self) {
        self.data[self.cur_chunk][self.cur_entry] = None;
        self.advance();
    }

    pub fn advance(&mut self) {
        self.cur_entry += 1;

        // have enough items to write a batch.
        if self.cur_entry >= self.data[self.cur_chunk].len() {
            self.yield_data_chunk();
        }

        // we wrote the last batch of the chunk, so make a new one.
        if self.cur_chunk >= self.data.len() {
            self.alloc_chunk();
        }
    }

    pub fn yield_data_chunk(&mut self) {
        let batch = std::mem::take(&mut self.data[self.cur_chunk]);

        for mut item in batch.into_iter().flatten() {
            let _ = item.write(&mut self.outbuf);
        }

        self.output.write(&self.outbuf);
        self.outbuf.clear();

        self.cur_chunk += 1;
        self.cur_entry = 0;
    }

    pub fn flush(&mut self) {
        let batch = std::mem::take(&mut self.data[self.cur_chunk]);

        for (i, entry) in batch.into_iter().enumerate() {
            if i >= self.cur_entry {
                break;
            }

            if let Some(mut dat) = entry {
                let _ = dat.write(&mut self.outbuf);
            }
        }

        self.output.write(&self.outbuf);
        self.outbuf.clear();

        self.cur_chunk += 1;
        self.cur_entry = 0;
    }
}

/// Defines how to get output data from iterators from a thread. If using a single thread, we can just print directly and not waste memory queueing output.
pub enum OutputMethod<W: std::io::Write, T: OrderedPileupOutput> {
    WriteDirectly(W),
    QueueForOutput(PileupOutputArray<T>),
}
