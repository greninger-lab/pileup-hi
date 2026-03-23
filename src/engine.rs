use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{IntervalJob, IntervalJobs, OrderedPileupOutput, OutputMethod, PileupOutputArray},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    utils::{get_writer_multi, OutputWriter},
};

use log::debug;
use std::sync::{Arc, Condvar, Mutex};

const OUTPUT_ARRAY_YIELD_SIZE: i64 = 2000;
pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

use anyhow::Error;
use log::{info, warn};
use std::io::BufWriter;

pub struct ThreadSignal {
    lock: Mutex<usize>,
    cvar: Condvar,
}

impl ThreadSignal {
    pub fn wait_while(&self) {
        std::mem::drop(
            self.cvar
                .wait_while(self.lock.lock().unwrap(), |free| *free == 0)
                .unwrap(),
        )
    }

    pub fn mark_running(&self) {
        *self.lock.lock().unwrap() -= 1;
        self.cvar.notify_one();
    }

    pub fn mark_done(&self) {
        *self.lock.lock().unwrap() += 1;
        self.cvar.notify_one();
    }
}

pub struct PileupWorker {
    jobid: usize,
    handle: Option<std::thread::JoinHandle<()>>,
    notify: Arc<ThreadSignal>,
}

impl PileupWorker {
    pub fn new(notify: Arc<ThreadSignal>) -> Self {
        Self {
            jobid: 0,
            handle: None,
            notify,
        }
    }

    fn is_finished(&mut self) -> bool {
        if let Some(ref handle) = self.handle {
            if handle.is_finished() {
                self.handle.take().unwrap().join().unwrap();
                return true;
            } else {
                return false;
            }
        }
        true
    }

    pub fn run<T>(
        &mut self,
        id: usize,
        params: PileupParams,
        job: IntervalJob,
        src: BamDataSource,
        o: T,
        out: OutputWriter,
        refseq: RefSeqHandle,
    ) where
        T: OrderedPileupOutput + 'static,
    {
        self.jobid = id;
        let notify = Arc::clone(&self.notify);

        self.handle = Some(std::thread::spawn(move || {
            notify.mark_running();

            let mut iterator = PileupIterator::new(
                &src,
                refseq,
                &params,
                o,
                OutputMethod::QueueForOutput(PileupOutputArray::new(
                    std::cmp::min((job.interval.len() / 10).max(1), OUTPUT_ARRAY_YIELD_SIZE) as usize,
                    out,
                )),
            )
            .unwrap();

            iterator.auto_loop2(&job.interval).unwrap();

            // signal that we're done.
            *job.done.lock().unwrap() = true;
            notify.mark_done();
        }));
    }
}

pub struct ThreadPool {
    workers: Vec<PileupWorker>,
    notify: Arc<ThreadSignal>,
}

impl ThreadPool {
    pub fn new(n_threads: usize) -> Self {
        let notify = Arc::new(ThreadSignal {
            lock: Mutex::new(n_threads),
            cvar: Condvar::new(),
        });

        let mut s = Self {
            notify,
            workers: Vec::with_capacity(n_threads),
        };

        (0..n_threads).for_each(|_| s.workers.push(PileupWorker::new(Arc::clone(&s.notify))));

        s
    }

    pub fn get_available(&mut self) -> Option<&mut PileupWorker> {
        self.notify.wait_while();

        for worker in self.workers.iter_mut() {
            if worker.is_finished() {
                return Some(worker);
            }
        }

        None
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
    refseq: RefSeq,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(in_params: InputParams, plp_params: PileupParams, output: T) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;
        let dest = OutputDataDest::from_string(&plp_params.output);

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = in_params.region {
            create_region_queue(&region, header)?
        } else {
            intervals_from_header(header)?
        };

        let refseq = RefSeq::new();

        Ok(Self {
            intervals,
            plp_params,
            src,
            output,
            dest,
            refseq,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.intervals.is_empty() {
            return Ok(());
        }

        // remove old output file if it exists.
        if let OutputDataDest::File(ref f) = self.dest {
            if std::fs::exists(f)? {
                warn!("Output file {} already exists! Overwriting...", f);

                if let Err(e) = std::fs::remove_file(f) {
                    warn!("Failed to remove file {f}; {e}. Output will be appended...");
                };
            }
        }

        if self.src.has_index()? {
            info!("Found index for for input file {}", self.src.fname()?);
        }

        if self.plp_params.threads == 1 {
            self.run_single()
        } else if !self.src.has_index()? {
            warn!(
                "User asked for more than {} threads but file is unindexed. Running in single-thread mode...",
                self.plp_params.threads
            );
            self.run_single()
        } else {
            info!("Running with {} threads...", self.plp_params.threads);
            self.run_multi()
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals.iter() {
            let main_writer: Box<dyn std::io::Write> = match self.dest {
                OutputDataDest::File(_) => Box::new(get_writer_multi(&self.dest, BUFWRITER_CAP, true, false)?),
                OutputDataDest::Stdout => Box::new(BufWriter::with_capacity(BUFWRITER_CAP, std::io::stdout().lock())),
            };

            let refseq_handle = self
                .refseq
                .yield_handle(&interval.name, self.plp_params.refseq.as_deref())?;

            let mut iterator = PileupIterator::new(
                &self.src,
                refseq_handle,
                &self.plp_params,
                self.output.clone(),
                OutputMethod::WriteDirectly(self.output.clone(), main_writer),
            )?;

            iterator.auto_loop2(interval)?;
        }
        Ok(())
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    /// The problem with this: all threads block until the last per-ref thread finishes.
    /// We need this to dynamically determinate how many chunks across ALL threads, dole out chunks
    /// in order to threads, and make sure file manager respects order.
    ///
    /// Thoughts: we need to be finish as many parts of the same reference before moving on to the
    /// next one. if we are processing eukaryotic genomes and have multiple chromosomes in memory, that
    /// will get ugly very fast; better to have one or two in memory at a time.
    pub fn run_multi(self) -> Result<(), Error> {
        let main_writer = get_writer_multi(&self.dest, BUFWRITER_CAP, true, false)?;

        let mut jobs = IntervalJobs::new(
            &self.intervals,
            self.plp_params.coords_per_thread,
            self.plp_params.threads as i64,
            main_writer,
        );

        let mut pool = ThreadPool::new(self.plp_params.threads);
        let mut n_jobs = 0;

        while !jobs.is_completed() {
            jobs.merge_completed()?;

            if let Some(worker) = pool.get_available() {
                if let Some(job) = jobs.queue.pop_front() {
                    n_jobs += 1;

                    let refseq_handle = self
                        .refseq
                        .yield_handle(&job.interval.name, self.plp_params.refseq.as_deref())?;

                    let writer = get_writer_multi(&job.out, BUFWRITER_CAP, true, false)?;

                    worker.run(
                        n_jobs,
                        self.plp_params.clone(),
                        job,
                        self.src.clone(),
                        self.output.clone(),
                        writer,
                        refseq_handle,
                    );
                }
            }
        }
        jobs.conclude()
    }
}
