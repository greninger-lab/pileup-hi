use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    errors::{Error, ErrorKind},
    jobqueue::IntervalJobs,
    output::{OrderedPileupOutput, OutputDestination, OutputFormat},
    params::{InputParams, PileupParams},
    pileup_iterator::{PileupIterator, PileupIteratorCore},
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    threading::ThreadPool,
    utils::get_writer_multi,
};

use log::{info, warn};
use std::sync::Arc;
use std::{
    cell::{Ref, RefCell},
    marker::PhantomData,
};

pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

struct PileupEngineQuery {
    intervals: Vec<GenomeInterval>,
    src: BamDataSource,
}

impl TryFrom<InputParams> for PileupEngineQuery {
    type Error = Error;
    fn try_from(value: InputParams) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&value.file)?;
        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = value.region {
            create_region_queue(&region, header)?
        } else {
            intervals_from_header(header)?
        };

        Ok(Self { intervals, src })
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    query: Option<RefCell<PileupEngineQuery>>,
    plp_params: PileupParams,
    _t: PhantomData<T>,
    dest: Option<OutputDataDest>,
    threads: usize,
    refseq: Option<RefSeq>,
}

/// An interface to generate pileups in memory, not writing to file. Single-threaded; parallelism is left up to the user.
pub struct PileupStream<T: OrderedPileupOutput + 'static> {
    engine: PileupEngine<T>,
}

impl<T: OrderedPileupOutput + 'static> PileupStream<T> {
    /// Return an iterator of pileups across the coordinates of the input regions specified. One iterator will be returned per region.
    pub fn get_iter(&mut self, input: InputParams) -> Result<Vec<PileupIterator<T>>, Error> {
        self.engine.submit(input)?;
        self.engine.yield_iterator()
    }
}

/// A pileup engine used to emit to files. Mulithreaded.
pub struct PileupSink<T: OrderedPileupOutput + 'static> {
    engine: PileupEngine<T>,
}

impl<T: OrderedPileupOutput + 'static> PileupSink<T> {
    /// Tell the engine to run over the specified input region
    pub fn submit(&mut self, input: InputParams) -> Result<(), Error> {
        self.engine.submit(input)
    }

    /// Run the engine (this function blocks)
    pub fn run(&self) -> Result<(), Error> {
        self.engine._run()
    }
}

///////////////////////////////////////////////////////////////

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    fn submit(&mut self, input: InputParams) -> Result<(), Error> {
        self.query = Some(RefCell::new(input.try_into()?));
        Ok(())
    }

    fn get_query(&self) -> Option<Ref<'_, PileupEngineQuery>> {
        self.query.as_ref().map(|b| b.borrow())
    }

    fn yield_iterator(&self) -> Result<Vec<PileupIterator<T>>, Error> {
        if let Some(query) = self.get_query() {
            let mut ret = Vec::with_capacity(query.intervals.len());

            for interval in query.intervals.iter() {
                let mut _iterator = PileupIteratorCore::new(
                    &query.src,
                    self.get_refseq(&query.intervals[0].name)?,
                    &self.plp_params,
                    OutputFormat::new(T::new(), OutputDestination::Memory),
                )?;

                _iterator.set_ref(interval.clone())?;

                ret.push(PileupIterator::from_iterator(_iterator))
            }

            Ok(ret)
        } else {
            Err(Error::from(ErrorKind::BadInputRegions(
                "user asked for iterator but no regions were loaded with submit()".to_string(),
            )))
        }
    }

    fn init_core(plp_params: PileupParams) -> Result<PileupEngine<T>, Error> {
        let refseq = if let Some(ref fasta) = plp_params.refseq {
            if !std::fs::exists(std::path::Path::new(fasta))? {
                return Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "Fasta file {fasta} doesn't exist!",
                ))));
            }
            Some(RefSeq::new(fasta.clone()))
        } else {
            None
        };

        Ok(Self {
            query: None,
            plp_params,
            _t: PhantomData,
            dest: None,
            threads: 1,
            refseq,
        })
    }

    /// Return an engine for iterating over records via the advance() API.
    pub fn init_stream(plp_params: PileupParams) -> Result<PileupStream<T>, Error> {
        Ok(PileupStream {
            engine: Self::init_core(plp_params)?,
        })
    }

    /// Return an engine for writing to FILE, as opposed to memory. The number of threads dictates the number of regions processed in parallel.
    pub fn init_sink(plp_params: PileupParams, output: &str, threads: usize) -> Result<PileupSink<T>, Error> {
        assert!(threads > 0, "invalid number of threads passed: {}", threads);
        let mut engine = Self::init_core(plp_params)?;
        engine.threads = threads;
        engine.dest = Some(OutputDataDest::from_string(output));
        Ok(PileupSink { engine })
    }

    fn _run(&self) -> Result<(), Error> {
        if let Some(query) = self.get_query() {
            // remove old output file if it exists.
            if let OutputDataDest::File(ref f) = self.dest.as_ref().unwrap() {
                if std::fs::exists(f)? {
                    warn!("Output file {} already exists! Overwriting...", f);

                    if let Err(e) = std::fs::remove_file(f) {
                        warn!("Failed to remove file {f}; {e}. Output will be appended...");
                    };
                }
            }

            if query.src.has_index()? {
                info!("Found index for for input file {}", query.src.fname()?);
            }

            if self.threads == 1 {
                self.run_all_1t()?;
            } else if !query.src.has_index()? {
                warn!(
                    "User asked for more than {} threads but file is unindexed. Running in single-thread mode...",
                    self.threads
                );
                self.run_all_1t()?;
            } else {
                info!("Running with {} threads...", self.threads);
                self.run_all_par()?;
            }
        };

        Ok(())
    }

    fn get_refseq(&self, ref_name: &str) -> Result<RefSeqHandle, Error> {
        if let Some(ref refseq) = self.refseq {
            refseq.yield_handle(ref_name)
        } else {
            Ok(Arc::new(None))
        }
    }

    /// Use a single thread for both processing and writing.
    fn run_all_1t(&self) -> Result<(), Error> {
        if let Some(ref query) = self.get_query() {
            for interval in query.intervals.iter() {
                let main_writer = get_writer_multi(self.dest.as_ref().unwrap(), BUFWRITER_CAP, true, false)?;

                let refseq_handle = self.get_refseq(&interval.name)?;

                let mut iterator = PileupIteratorCore::new(
                    &query.src,
                    refseq_handle,
                    &self.plp_params,
                    OutputFormat::new(T::new(), OutputDestination::Writer(main_writer)),
                )?;

                iterator.auto_loop2(interval)?;
            }
        }
        Ok(())
    }

    /// Split up a list of input genomic intervals into smaller chunks to be processed in parallel. Chunks are first written to temporary output files before being merged into the user-specified output file.
    fn run_all_par(&self) -> Result<(), Error> {
        if let Some(query) = self.get_query() {
            let mut jobs = IntervalJobs::new(
                &query.intervals,
                self.plp_params.coords_per_thread,
                self.threads as i64,
                self.dest.as_ref().unwrap().clone(),
            );

            let mut pool = ThreadPool::new(self.threads);
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
                            query.src.clone(),
                            T::new(),
                            refseq_handle,
                        );
                    }
                }
            }
            jobs.conclude()?;
        };
        Ok(())
    }
}
