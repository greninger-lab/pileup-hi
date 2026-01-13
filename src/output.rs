use crate::alignment::PileupAlignment;
use crate::bamio::OutputDataDest;
use crate::engine::BUFWRITER_CAP;
use crate::utils::{get_writer, temp_fname, OutputWriter};
use anyhow::Error;
use log::warn;
use std::fs::File;
use std::io::{BufReader, Write};
use std::sync::Mutex;

pub static FILE_MERGE_SINGLETON: Mutex<OutputFileMerge> = Mutex::new(OutputFileMerge {
    outfile: OutputDataDest::Stdout,
    subfiles: vec![],
});

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Get the reference of the pileup
    #[allow(dead_code)]
    fn tid(&self) -> i32;
    /// Get the coordinate of the pileup
    #[allow(dead_code)]
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

/// Used to keep track of our main output file and the subfiles we want to merge into it. Subfiles
/// are ordered by thread ID.
#[derive(Clone)]
pub struct OutputFileMerge {
    pub outfile: OutputDataDest,
    pub subfiles: Vec<OutputDataDest>,
}

impl OutputFileMerge {
    /// Copy the data in the subfiles over to the main file
    pub fn merge<W: std::io::Write>(&self, mut dest: W) -> Result<(), Error> {
        for s in &self.subfiles {
            match s {
                OutputDataDest::Stdout => anyhow::bail!("cannot merge from stdout! Critical error"),
                OutputDataDest::File(ref f) => {
                    let fhandle = File::open(f)?;
                    let mut reader = BufReader::with_capacity(2 * 1024 * 1024, fhandle);
                    std::io::copy(&mut reader, &mut dest)?;
                    if let Err(e) = std::fs::remove_file(f) {
                        match e.kind() {
                            std::io::ErrorKind::NotFound => (),
                            _ => anyhow::bail!(e),
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// Delete all files at once, useful if we abruptly exit
    pub fn cleanup(&mut self) -> Result<(), Error> {
        for s in &self.subfiles {
            if let OutputDataDest::File(f) = s {
                if let Err(e) = std::fs::remove_file(f) {
                    match e.kind() {
                        std::io::ErrorKind::NotFound => (),
                        _ => anyhow::bail!(e),
                    }
                }
            }
        }

        Ok(())
    }

    /// If main thread, return the writer to the final output file
    pub fn get_writer(&self, thread_idx: usize) -> Result<OutputWriter, Error> {
        if thread_idx == 0 {
            get_writer(&self.outfile, BUFWRITER_CAP, true, true)
        } else {
            get_writer(&self.subfiles[thread_idx - 1], BUFWRITER_CAP, true, false)
        }
    }
}

pub fn generate_subfile_dests(outprefix: &str, n: usize, suffix: &str) -> Vec<OutputDataDest> {
    let mut ret = Vec::with_capacity(n);
    for i in 0..n {
        let temp = temp_fname(outprefix, &i.to_string(), suffix);
        ret.push(OutputDataDest::File(temp));
    }

    ret
}

pub fn setup_exit_handler() {
    ctrlc::set_handler(|| {
        warn!("Received termination signal. Cleaning up intermediate files...");
        if let Ok(mut outputs) = FILE_MERGE_SINGLETON.lock() {
            outputs
                .cleanup()
                .expect("Failed to cleanup temp files during termination sequence");
        }

        std::process::exit(130);
    })
    .expect("Failed to set exit handler")
}

pub struct PileupOutputArray<T: OrderedPileupOutput> {
    data: Vec<T>,
    writable: Vec<bool>,
    cur: usize,
    capacity: usize,
    writer: OutputWriter,
}

impl<T: OrderedPileupOutput> PileupOutputArray<T> {
    pub fn new(capacity: usize, writer: OutputWriter) -> Self {
        Self {
            data: vec![T::new(); capacity],
            writable: vec![true; capacity],
            cur: 0,
            capacity,
            writer,
        }
    }

    pub fn cur_mut(&mut self) -> &mut T {
        &mut self.data[self.cur]
    }

    // no-op
    pub fn push(&mut self) {}

    pub fn tombstone(&mut self) {
        self.writable[self.cur] = false
    }

    pub fn advance(&mut self) -> Result<(), Error> {
        self.cur += 1;

        if self.cur >= self.capacity {
            self.write_all()?;
        }

        Ok(())
    }

    pub fn write_all(&mut self) -> Result<(), Error> {
        for i in 0..self.cur {
            if self.writable[i] {
                self.data[i].write(&mut self.writer)?;
            } else {
                self.data[i].clear();
            }
        }

        self.cur = 0;
        self.writable.fill(true);
        Ok(())
    }
}

/// Defines how to get output data from iterators from a thread. If using a single thread, we can just print directly and not waste memory queueing output.
pub enum OutputMethod<T: OrderedPileupOutput> {
    WriteDirectly(Box<dyn Write>),
    QueueForOutput(PileupOutputArray<T>),
}
