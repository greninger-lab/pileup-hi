use crate::alignment::PileupAlignment;
use crate::bamio::OutputDataDest;
use crate::engine::BUFWRITER_CAP;
use crate::utils::{get_writer_multi, temp_fname, OutputWriter};
use crate::{position_queue::GenomeInterval, refseq::RefSeqHandle};
use anyhow::Error;
use crossbeam::channel::{unbounded, Receiver, Sender};
use log::{info, warn};
use std::collections::VecDeque;
use std::fs::File;
use std::io::BufReader;
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

    #[allow(dead_code)]
    fn new() -> Self;
}

/// A job a worker thread can accept. Contains temp output file handle, the interval to process, and a "done" state.
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

/// Data structure storing tasks to send to threads. Each interval job (chunk) is mapped to the larger interval it was fragmented from. Once all chunks are completed, their respective temp files are copied into the main output file and deleted. Importantly, larger intervals are ORDERED.
pub struct IntervalJobs {
    map: VecDeque<(GenomeInterval, Vec<IntervalJob>)>,
    pub queue: VecDeque<IntervalJob>,
    handle: std::thread::JoinHandle<()>,
    s: Sender<Vec<IntervalJob>>,
}

impl IntervalJobs {
    pub fn new(intervals: &[GenomeInterval], min_coords_per_thread: i64, threads: i64, output: OutputDataDest) -> Self {
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
            let mut main_writer = get_writer_multi(&output, BUFWRITER_CAP, true, false).unwrap();
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

    /// Check if we have any intervals with all chunks completed.
    pub fn merge_completed(&mut self) -> Result<(), Error> {
        let mut done = 0;

        if let Some((interval, pending)) = self.map.front() {
            for tmp in pending {
                if *tmp.done.lock().unwrap() {
                    done += 1;
                }
            }

            assert!(done <= pending.len());

            // all chunks have been marked "done" by their assigned workers
            if done == pending.len() {
                info!("Finished ref {}", interval.name);
                let (_, to_merge) = self.map.pop_front().unwrap();
                self.s.send(to_merge)?;
            }
        }

        Ok(())
    }

    /// Should be called when we've done everything (e.g. self.is_completed()). Does one final merge
    /// and signals the writer thread to finish.
    pub fn conclude(mut self) -> Result<(), Error> {
        self.merge_completed()?;
        drop(self.s);
        self.handle.join().expect("Failed to join writer thread");
        Ok(())
    }
}

/// Tell the program on unexpected exit (e.g. SIGTERM, ctrl-c) to delete all temp files it created and hasn't yet merged to final output
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
pub struct OutputFormat<T: OrderedPileupOutput> {
    output: T,
    writer: OutputWriter,
}

impl<T: OrderedPileupOutput> OutputFormat<T> {
    pub fn new(output: T, writer: OutputWriter) -> Self {
        Self { output, writer }
    }

    pub fn reject(&mut self) -> bool {
        self.output.clear();
        false
    }

    pub fn cur(&mut self) -> &mut T {
        &mut self.output
    }

    pub fn take(&mut self) -> Result<bool, Error> {
        self.output.write(&mut self.writer)?;
        Ok(true)
    }

    pub fn check(&mut self, emit: bool) -> Result<bool, Error> {
        if emit {
            self.take()
        } else {
            Ok(self.reject())
        }
    }
}
