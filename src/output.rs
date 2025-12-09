#![allow(dead_code)]
use crate::alignment::PileupAlignment;
use anyhow::Error;
use crossbeam::channel::{bounded, Receiver, Sender};
use std::io::BufWriter;
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

/// Defines how to get output data from iterators from a thread. If using a single thread, we can just print directly and not waste memory queueing output.
/// have to care about queue-ing output.
pub enum OutputMethod<W: std::io::Write, T: OrderedPileupOutput> {
    WriteDirectly(W),
    QueueForOutput(Vec<T>),
}

////////////////
// Begin defs for PileupOutputAggregator
////////////////

pub struct PileupOutputAggregator<T>
where
    T: OrderedPileupOutput,
{
    pub input_handle: Option<Sender<T>>,
    pub join_handle: Option<JoinHandle<()>>,
}

impl<T: OrderedPileupOutput + 'static> PileupOutputAggregator<T> {
    pub fn new() -> Self {
        Self {
            input_handle: None,
            join_handle: None,
        }
    }

    pub fn get_output_handle(&self) -> Option<Sender<T>> {
        self.input_handle.clone()
    }

    pub fn terminate(self) -> Result<(), Error> {
        match (self.input_handle, self.join_handle) {
            (None, _) | (_, None) => anyhow::bail!("attempted to terminate an unitialized aggregator."),
            (Some(snd), Some(join_handle)) => {
                drop(snd);
                join_handle.join().expect("failed to join output aggregator");
                Ok(())
            }
        }
    }

    pub fn run(&mut self) {
        let (s, r): (Sender<T>, Receiver<T>) = bounded(10_000_000);
        let j = std::thread::spawn(move || {
            let mut writer = BufWriter::new(std::io::stdout().lock());
            r.into_iter().for_each(|mut o| o.write(&mut writer).unwrap());
        });

        self.join_handle = Some(j);
        self.input_handle = Some(s);
    }
}
