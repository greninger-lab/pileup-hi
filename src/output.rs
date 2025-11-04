use crate::threading::PileupWorkerState;
use anyhow::Error;
use crossbeam::channel::{unbounded, Receiver, Sender};
const PILEUP_OUTPUT_BUF_PURGE_THRES: usize = 3;

pub trait OrderedPileupOutput {
    fn tid(&self) -> i32;
    fn pos(&self) -> i64;
    fn write(&mut self) -> Result<(), Error>;
}

pub enum PileupOutputState<T: OrderedPileupOutput> {
    Closed,
    Open(Sender<T>),
}

pub struct PileupOutputAggregator<T>
where
    T: OrderedPileupOutput,
{
    pub input_state: PileupOutputState<T>,
    pub worker_state: PileupWorkerState,
}

impl<T: OrderedPileupOutput + Send + 'static> PileupOutputAggregator<T> {
    pub fn new() -> Self {
        Self {
            input_state: PileupOutputState::Closed,
            worker_state: PileupWorkerState::Off,
        }
    }

    pub fn get_output_handle(&self) -> Option<Sender<T>> {
        match &self.input_state {
            PileupOutputState::Closed => None,
            PileupOutputState::Open(s) => Some(s.clone()),
        }
    }

    pub fn terminate(self) -> Result<(), Error> {
        match (self.input_state, self.worker_state) {
            (PileupOutputState::Closed, _) | (_, PileupWorkerState::Off) => {
                anyhow::bail!("Cannot terminate an output channel that never started!")
            }

            (PileupOutputState::Open(s), PileupWorkerState::Running(j)) => {
                drop(s);
                j.join().unwrap();
                Ok(())
            }
        }
    }

    pub fn run(&mut self) {
        let (s, r): (Sender<T>, Receiver<T>) = unbounded();
        let j = std::thread::spawn(move || {
            let mut output_queue: Vec<T> = Vec::with_capacity(PILEUP_OUTPUT_BUF_PURGE_THRES);
            let mut next_expected_order = 0;

            // initialize next expected order
            if let Ok(mut out) = r.recv() {
                next_expected_order = out.pos();
                out.write().unwrap();
            }

            while let Ok(mut out) = r.recv() {
                if out.pos() == next_expected_order + 1 {
                    next_expected_order += 1;
                    out.write().unwrap();
                    continue;
                }

                output_queue.push(out);

                if output_queue.len() >= PILEUP_OUTPUT_BUF_PURGE_THRES {
                    output_queue
                        .sort_by(|a, b| a.tid().cmp(&b.tid()).then_with(|| a.pos().cmp(&b.pos())));

                    let mut processable_count = 0;

                    for item in &output_queue {
                        if item.pos() == next_expected_order + 1 {
                            processable_count += 1;
                            next_expected_order += 1;
                        }
                    }

                    for mut out in output_queue.drain(..processable_count) {
                        out.write().unwrap();
                    }

                    // output_queue.shrink_to(0);
                }
            }
        });

        self.worker_state = PileupWorkerState::Running(j);
        self.input_state = PileupOutputState::Open(s.clone());
    }
}
