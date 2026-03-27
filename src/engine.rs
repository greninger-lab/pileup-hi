use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{IntervalJob, IntervalJobs, OrderedPileupOutput, OutputFormat},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    utils::get_writer_multi,
};

use std::sync::{Arc, Condvar, Mutex};

use anyhow::Error;
use log::{info, warn};

pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

////////////
// THREADING
////////////

// A conditional variable that tracks the number of threads available to accept a genomic interval
// job.
pub struct ThreadSignal {
    lock: Mutex<usize>, // number of available threads
    cvar: Condvar,      // what to ping when number of free threads changes
}

impl ThreadSignal {
    // Wait until there is at least one available thread
    pub fn wait_while(&self) {
        std::mem::drop(
            self.cvar
                .wait_while(self.lock.lock().unwrap(), |free| *free == 0)
                .unwrap(),
        )
    }

    // Notify that a thread has started working and become unavailable
    pub fn mark_running(&self) {
        *self.lock.lock().unwrap() -= 1;
        self.cvar.notify_one();
    }

    // Notify that a thread has finished its job and is now available
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
        refseq: RefSeqHandle,
    ) where
        T: OrderedPileupOutput + 'static,
    {
        self.jobid = id;
        let notify = Arc::clone(&self.notify);

        self.handle = Some(std::thread::spawn(move || {
            notify.mark_running();

            let out = get_writer_multi(&job.out, BUFWRITER_CAP, true, false).unwrap();

            let mut iterator = PileupIterator::new(&src, refseq, &params, OutputFormat::new(o, out)).unwrap();

            iterator.auto_loop2(&job.interval).unwrap();

            // signal that we're done.
            *job.done.lock().unwrap() = true;
            notify.mark_done();
        }));
    }
}

/// Very simple abstraction over a collection of threads that can take jobs.
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

        self.workers.iter_mut().find_map(|w| w.is_finished().then_some(w))
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
    refseq: Option<RefSeq>,
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

        let refseq = if let Some(ref fasta) = plp_params.refseq {
            if !std::fs::exists(std::path::Path::new(fasta))? {
                anyhow::bail!("Fasta file provided ({}), doesn't exist!", fasta);
            }
            Some(RefSeq::new(fasta.clone()))
        } else {
            None
        };

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

    fn get_refseq(&self, ref_name: &str) -> Result<RefSeqHandle, Error> {
        if let Some(ref refseq) = self.refseq {
            refseq.yield_handle(ref_name)
        } else {
            Ok(Arc::new(None))
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals.iter() {
            let main_writer = get_writer_multi(&self.dest, BUFWRITER_CAP, true, false)?;

            let refseq_handle = self.get_refseq(&interval.name)?;

            let mut iterator = PileupIterator::new(
                &self.src,
                refseq_handle,
                &self.plp_params,
                OutputFormat::new(self.output.clone(), main_writer),
            )?;

            iterator.auto_loop2(interval)?;
        }
        Ok(())
    }

    /// Split up a list of input genomic intervals into smaller chunks to be processed in parallel. Chunks are first written to temporary output files before being merged into the user-specified output file.
    pub fn run_multi(self) -> Result<(), Error> {
        let mut jobs = IntervalJobs::new(
            &self.intervals,
            self.plp_params.coords_per_thread,
            self.plp_params.threads as i64,
            self.dest.clone(),
        );

        let mut pool = ThreadPool::new(self.plp_params.threads);
        let mut n_jobs = 0;

        while !jobs.is_completed() {
            jobs.merge_completed()?;

            if let Some(worker) = pool.get_available() {
                if let Some(job) = jobs.queue.pop_front() {
                    n_jobs += 1;

                    let refseq_handle = self.get_refseq(&job.interval.name)?;

                    worker.run(
                        n_jobs,
                        self.plp_params.clone(),
                        job,
                        self.src.clone(),
                        self.output.clone(),
                        refseq_handle,
                    );
                }
            }
        }
        jobs.conclude()
    }
}
