use crate::{
    bamio::{BamDataSource, BamReader},
    output::{OrderedPileupOutput, OutputMethod, PileupOutputAggregator},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, GenomeInterval, PositionQueue},
};

use anyhow::Error;
use rayon::prelude::*;
use std::{collections::VecDeque, io::BufWriter};

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
        Self { interval, params, src }
    }

    pub fn run<T>(&mut self, o: T) -> Vec<T>
    where
        T: OrderedPileupOutput + 'static,
    {
        let iterator = PileupIterator::new(
            &self.src,
            &self.params,
            o,
            OutputMethod::<DummyOutputWriter, T>::QueueForOutput(Vec::with_capacity(10_000)),
        )
        .unwrap();

        iterator
            ._auto_loop_yield_batch(&PositionQueue {
                queue: vec![self.interval.clone()],
            })
            .unwrap()
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: PositionQueue,
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
            PositionQueue::new(header)?
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
        let lock = BufWriter::new(std::io::stdout().lock());
        let mut iterator = PileupIterator::new(
            &self.src,
            &self.plp_params,
            self.output,
            OutputMethod::WriteDirectly(lock),
        )?;
        iterator._auto_loop_output_each(&self.intervals)
    }

    /// Use separate threads for processing and writing. Each processing thread owns its own reads into the BAM file, index, and any other files.
    pub fn run_multi(self) -> Result<(), Error> {
        for interval in self.intervals.queue {
            let mut agg: PileupOutputAggregator<T> = PileupOutputAggregator::new();
            agg.run();
            let output_handle = agg.get_output_handle().unwrap();

            let subintervals = interval.chunks(1_000_000).collect::<VecDeque<GenomeInterval>>();

            let threadpool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.threads)
                .build()
                .unwrap();
            let src = &self.src.clone();

            // thank you Seth Stadick for this this blazingly-fast rayon usage pattern.
            threadpool.install(|| {
                subintervals
                    .par_iter()
                    .flat_map(|chunk| {
                        let mut worker = PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
                        worker.run(self.output.clone())
                    })
                    .for_each(|o| {
                        output_handle.send(o).unwrap();
                    });
            });

            drop(output_handle);
            agg.terminate()?;
        }

        Ok(())
    }
}
