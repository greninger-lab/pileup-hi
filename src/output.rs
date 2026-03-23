use crate::alignment::PileupAlignment;
use crate::bamio::OutputDataDest;
use crate::utils::{temp_fname, OutputWriter};
use crate::{position_queue::GenomeInterval, refseq::RefSeqHandle};
use anyhow::Error;
use crossbeam::channel::{unbounded, Receiver, Sender};
use log::{info, warn};
use std::collections::VecDeque;
use std::fs::File;
use std::io::{BufReader, Write};
use std::sync::{Arc, Mutex};

pub static FILE_MERGE_SINGLETON: Mutex<Vec<OutputDataDest>> = Mutex::new(vec![]);

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
    fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error>;
    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, refseq: &RefSeqHandle);
    fn write<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), Error>;
    fn depth(&self) -> u32;
    fn clear(&mut self);
    fn new() -> Self;
}

pub struct IntervalJobInner {
    pub out: OutputDataDest,
    pub interval: GenomeInterval,
    pub done: Mutex<bool>,
}

impl IntervalJobInner {
    fn new(interval: &GenomeInterval) -> Self {
        Self {
            out: OutputDataDest::from_string(&temp_fname(
                &format!("{}:{}-{}", interval.name, interval.start, interval.end),
                "",
                ".temp",
            )),
            done: Mutex::new(false),
            interval: interval.clone(),
        }
    }
}

pub type IntervalJob = Arc<IntervalJobInner>;

pub struct IntervalJobs {
    map: VecDeque<(GenomeInterval, Vec<IntervalJob>)>,
    pub queue: VecDeque<IntervalJob>,
    handle: std::thread::JoinHandle<()>,
    s: Sender<Vec<IntervalJob>>,
}

impl IntervalJobs {
    pub fn new(
        intervals: &[GenomeInterval],
        min_coords_per_thread: i64,
        threads: i64,
        mut main_writer: OutputWriter,
    ) -> Self {
        let mut map: VecDeque<(GenomeInterval, Vec<IntervalJob>)> = VecDeque::new();
        let mut queue: VecDeque<IntervalJob> = VecDeque::new();
        let mut lock = FILE_MERGE_SINGLETON.lock().unwrap();

        for interval in intervals {
            let chunks = if interval.len() < min_coords_per_thread {
                interval.chunks(min_coords_per_thread)
            } else {
                interval.n_chunks(threads)
            }
            .map(|c| Arc::new(IntervalJobInner::new(&c)))
            .collect::<Vec<IntervalJob>>();

            chunks.iter().for_each(|c| {
                queue.push_back(c.clone());
                lock.push(c.out.clone());
            });

            map.push_back((interval.clone(), chunks.clone()));
        }

        let (s, r): (Sender<Vec<IntervalJob>>, Receiver<Vec<IntervalJob>>) = unbounded();

        let handle = std::thread::spawn(move || {
            while let Ok(temps) = r.recv() {
                for tmp in temps {
                    match tmp.out {
                        OutputDataDest::Stdout => panic!("cannot merge from stdout! Critical error"),
                        OutputDataDest::File(ref f) => {
                            match File::open(f) {
                                Err(e) => {
                                    match e.kind() {
                                        std::io::ErrorKind::NotFound => (),
                                        _ => panic!("Failed to open output file for merging: {}", e),
                                    };
                                }

                                Ok(f) => {
                                    let mut reader = BufReader::with_capacity(2 * 1024 * 1024, f);
                                    std::io::copy(&mut reader, &mut main_writer).unwrap();
                                }
                            }
                            if let Err(e) = std::fs::remove_file(f) {
                                match e.kind() {
                                    std::io::ErrorKind::NotFound => (),
                                    _ => panic!("{}", e),
                                }
                            }
                        }
                    }
                }
            }
        });

        Self { map, handle, queue, s }
    }

    pub fn is_completed(&self) -> bool {
        self.map.is_empty()
    }

    pub fn merge_completed(&mut self) -> Result<(), Error> {
        let mut done = 0;

        if let Some((interval, pending)) = self.map.front() {
            for tmp in pending {
                if *tmp.done.lock().unwrap() {
                    done += 1;
                }
            }

            assert!(done <= pending.len());
            if done == pending.len() {
                info!("Finished ref {}", interval.name);
                let (_, to_merge) = self.map.pop_front().unwrap();
                self.s.send(to_merge)?;
            }
        }

        Ok(())
    }

    pub fn conclude(mut self) -> Result<(), Error> {
        self.merge_completed()?;
        drop(self.s);
        self.handle.join().expect("Failed to join writer thread");
        Ok(())
    }
}

pub fn setup_exit_handler() {
    ctrlc::set_handler(|| {
        warn!("Received termination signal. Cleaning up intermediate files...");
        if let Ok(outputs) = FILE_MERGE_SINGLETON.lock() {
            for t in outputs.iter() {
                match t {
                    OutputDataDest::Stdout => (),
                    OutputDataDest::File(ref f) => {
                        if let Err(e) = std::fs::remove_file(f) {
                            match e.kind() {
                                std::io::ErrorKind::NotFound => (),
                                _ => eprintln!("{e}"),
                            }
                        }
                    }
                }
            }
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
    WriteDirectly(T, Box<dyn Write>),
    QueueForOutput(PileupOutputArray<T>),
}

impl<T: OrderedPileupOutput> OutputMethod<T> {
    pub fn cur(&mut self) -> &mut T {
        match self {
            Self::WriteDirectly(output, _writer) => output,
            Self::QueueForOutput(output_arr) => output_arr.cur_mut(),
        }
    }

    pub fn reject(&mut self) -> Result<bool, Error> {
        match self {
            Self::WriteDirectly(output, _writer) => output.clear(),
            Self::QueueForOutput(output_arr) => {
                output_arr.tombstone();
                output_arr.advance()?;
            }
        }
        Ok(false)
    }

    pub fn check(&mut self, emit: bool) -> Result<bool, Error> {
        if emit {
            self.take()?;
            Ok(true)
        } else {
            self.reject()?;
            Ok(false)
        }
    }

    pub fn take(&mut self) -> Result<bool, Error> {
        match self {
            Self::WriteDirectly(output, writer) => output.write(writer)?,
            Self::QueueForOutput(output_arr) => {
                output_arr.push();
                output_arr.advance()?;
            }
        }

        Ok(true)
    }
}
