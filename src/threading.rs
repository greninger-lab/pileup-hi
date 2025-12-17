use crate::{
    bamio::{BamDataSource, BamReader},
    output::{
        OrderedPileupOutput, OutputMethod, PileupOutputAggregator, PileupOutputArray,
        PileupOutputChunk,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
};

const OUTPUT_ARRAY_YIELD_SIZE: usize = 2000;

use anyhow::Error;
use crossbeam::channel::Sender;
use rayon::prelude::*;
use std::io::BufWriter;

pub struct PileupWorker {
    interval: GenomeInterval,
    params: PileupParams,
    src: BamDataSource,
}

pub struct DummyOutputWriter {}
impl std::io::Write for DummyOutputWriter {
    fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
        Ok(0)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl PileupWorker {
    pub fn new(params: PileupParams, interval: GenomeInterval, src: BamDataSource) -> Self {
        Self {
            interval,
            params,
            src,
        }
    }

    pub fn run<T>(&mut self, o: T, snd: Sender<PileupOutputChunk>, id: u8)
    where
        T: OrderedPileupOutput + 'static,
    {
        let iterator = PileupIterator::new(
            &self.src,
            &self.params,
            o,
            OutputMethod::<DummyOutputWriter, T>::QueueForOutput(
                PileupOutputArray::new(self.interval.len(), OUTPUT_ARRAY_YIELD_SIZE, snd, id)
                    .unwrap(),
            ),
        )
        .unwrap();

        iterator
            ._auto_loop_yield_batch(std::slice::from_ref(&self.interval))
            .unwrap()
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    threads: usize,
    output: T,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(
        in_params: InputParams,
        plp_params: PileupParams,
        threads: usize,
        output: T,
    ) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = &in_params.region {
            create_region_queue(region, header)?
        } else {
            // PositionQueue::new(header)?
            intervals_from_header(header)?
        };

        Ok(Self {
            intervals,
            plp_params,
            src,
            threads,
            output,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.intervals.len() == 1 || self.threads == 1 || !self.src.has_index()? {
            self.run_single()
        } else {
            self.run_multi()
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals {
            let lock = BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout().lock());
            let mut iterator = PileupIterator::new(
                &self.src,
                &self.plp_params,
                self.output.clone(),
                OutputMethod::WriteDirectly(lock),
            )?;

            iterator._auto_loop_output_each(&[interval])?;
        }
        Ok(())
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    pub fn run_multi(self) -> Result<(), Error> {
        for interval in self.intervals {
            let mut agg: PileupOutputAggregator = PileupOutputAggregator::new();
            agg.run(self.src.fname().unwrap(), self.threads);
            let output_handle = agg.get_output_handle().unwrap();

            let subintervals = interval.chunks(1_000_000).collect::<Vec<GenomeInterval>>();

            let threadpool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .unwrap();

            let src = &self.src.clone();

            // Process ALL subintervals in parallel at once
            threadpool.install(|| {
                subintervals.par_iter().for_each(|chunk| {
                    let thread_id = rayon::current_thread_index().unwrap() as u8;
                    let mut worker =
                        PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
                    worker.run(self.output.clone(), output_handle.clone(), thread_id);
                });
            });

            drop(output_handle);
            agg.terminate()?;
        }

        Ok(())
    }
}
