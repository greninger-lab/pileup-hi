use crate::{
    bamio::BamDataSource,
    engine::BUFWRITER_CAP,
    jobqueue::IntervalJob,
    output::{OrderedPileupOutput, OutputFormat},
    params::PileupParams,
    pileup_iterator::PileupIterator,
    refseq::RefSeqHandle,
    utils::get_writer_multi,
};

use std::sync::{Arc, Condvar, Mutex};

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
