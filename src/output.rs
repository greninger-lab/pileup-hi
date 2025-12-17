#![allow(dead_code)]
use crate::alignment::PileupAlignment;
use crate::utils::temp_fname;
use anyhow::Error;
use crossbeam::channel::{bounded, unbounded, Receiver, Sender};
use std::fs::File;
use std::io::{stdout, BufRead, BufReader, BufWriter, Read, Write};
use std::thread::JoinHandle;

const PILEUP_OUTPUT_BUF_PURGE_THRES: usize = 1_000_000;

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone {
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
    writer: BufWriter<File>,
}

impl TempOutputHandle {
    pub fn new(prefix: &str, id: &str, writer_cap: usize) -> Result<Self, Error> {
        let fname = temp_fname(prefix, id, "temp");
        let file = File::create(&fname)?;

        Ok(Self {
            fname,
            writer: BufWriter::with_capacity(writer_cap, file),
        })
    }

    pub fn write(&mut self, data: &[u8]) {
        self.writer.write_all(data);
    }

    pub fn delete(&mut self) -> Result<(), Error> {
        std::fs::remove_file(std::path::Path::new(&self.fname)).map_err(Error::msg)
    }
}

fn merge_temp_outputs<W: std::io::Write>(
    outputs: &mut [TempOutputHandle],
    mut dest: W,
) -> Result<(), Error> {
    for temp in outputs {
        let mut buffer = [0u8; 8192];
        let mut reader = BufReader::with_capacity(2 * 1024 * 1024, File::open(&temp.fname)?);
        loop {
            let bytes_read = reader.read(&mut buffer)?;
            if bytes_read == 0 {
                break;
            }

            dest.write_all(&buffer)?;
        }

        temp.delete();
    }

    dest.flush()?;

    Ok(())
}

/// a pre-allocated array of pileup outputs used to buffer when processing genome intervals (intended for multithreading)
pub struct PileupOutputArray<T: OrderedPileupOutput> {
    data: Vec<Vec<Option<T>>>,
    cur_chunk: usize,
    cur_entry: usize,
    chunk_size: usize,
    output: Sender<PileupOutputChunk>,
    id: u8,
    outbuf: Vec<u8>,
}

impl<T: OrderedPileupOutput> PileupOutputArray<T> {
    pub fn new(
        capacity: usize,
        chunks: usize,
        output: Sender<PileupOutputChunk>,
        id: u8,
    ) -> Result<Self, Error> {
        let chunk_size = capacity / chunks;
        let remainder = capacity % chunks;

        let mut data = Vec::with_capacity(chunks);

        for _ in 0..chunks - 1 {
            data.push(vec![Some(T::new()); chunk_size]);
        }

        let final_size = if remainder != 0 {
            chunk_size + remainder
        } else {
            chunk_size
        };

        data.push(vec![Some(T::new()); final_size]);
        let outbuf = Vec::with_capacity(chunk_size * size_of::<T>());

        Ok(Self {
            data,
            cur_entry: 0,
            cur_chunk: 0,
            chunk_size,
            output,
            outbuf,
            id,
        })
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

        if self.cur_entry >= self.data[self.cur_chunk].len() {
            self.yield_data_chunk(
                // (self.chunk_size * self.cur_chunk.saturating_sub(1)) as u32
                //     + self.cur_entry as u32
                //     + self.pos_start as u32,
            );
        }
    }

    pub fn yield_data_chunk(&mut self) {
        let batch = std::mem::take(&mut self.data[self.cur_chunk]);

        for mut item in batch.into_iter().flatten() {
            item.write(&mut self.outbuf);
        }

        let out = std::mem::replace(
            &mut self.outbuf,
            Vec::with_capacity(self.chunk_size * size_of::<T>()),
        );

        self.output
            .send(PileupOutputChunk {
                data: out,
                id: self.id as usize,
            })
            .unwrap();

        self.cur_chunk += 1;
        self.cur_entry = 0;
    }
}

/// Defines how to get output data from iterators from a thread. If using a single thread, we can just print directly and not waste memory queueing output.
/// have to care about queue-ing output.
pub enum OutputMethod<W: std::io::Write, T: OrderedPileupOutput> {
    WriteDirectly(W),
    QueueForOutput(PileupOutputArray<T>),
}

////////////////
// Begin defs for PileupOutputAggregator
////////////////

pub struct PileupOutputAggregator {
    pub input_handle: Option<Sender<PileupOutputChunk>>,
    pub join_handle: Option<JoinHandle<()>>,
}

impl PileupOutputAggregator {
    pub fn new() -> Self {
        Self {
            input_handle: None,
            join_handle: None,
        }
    }

    pub fn get_output_handle(&self) -> Option<Sender<PileupOutputChunk>> {
        self.input_handle.clone()
    }

    pub fn terminate(self) -> Result<(), Error> {
        match (self.input_handle, self.join_handle) {
            (None, _) | (_, None) => {
                anyhow::bail!("attempted to terminate an unitialized aggregator.")
            }
            (Some(snd), Some(join_handle)) => {
                drop(snd);
                join_handle
                    .join()
                    .expect("failed to join output aggregator");
                Ok(())
            }
        }
    }

    pub fn run(&mut self, outprefix: String, threads: usize) {
        let (s, r): (Sender<PileupOutputChunk>, Receiver<PileupOutputChunk>) = bounded(10000);

        let j = std::thread::spawn(move || {
            let mut writers = Vec::with_capacity(threads - 1);

            for i in 0..threads - 1 {
                let temp =
                    TempOutputHandle::new(&outprefix, &i.to_string(), 2 * 1024 * 1024).unwrap();
                writers.push(temp);
            }

            let mut main_writer = BufWriter::with_capacity(2 * 1024 * 1024, stdout().lock());

            r.into_iter().for_each(|o| {
                if o.id == 0 {
                    main_writer.write_all(&o.data);
                    // o.data
                    //     .into_iter()
                    //     .flatten()
                    //     .for_each(|mut pos| pos.write(&mut main_writer).unwrap());
                } else {
                    let tempout = &mut writers[o.id - 1];
                    tempout.write(&o.data);
                    // o.data
                    //     .into_iter()
                    //     .flatten()
                    //     .for_each(|pos| tempout.write(pos));
                }
            });

            main_writer.flush().unwrap();

            merge_temp_outputs(&mut writers, main_writer).unwrap();
        });

        self.join_handle = Some(j);
        self.input_handle = Some(s);
    }
}
