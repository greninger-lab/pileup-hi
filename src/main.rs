use anyhow::Error;
use clap::Parser;

mod read_buf;
mod rpileup;

#[derive(Parser)]
pub struct Args {
    pub input: String,
}

fn main() -> Result<(), Error> {
    let args = Args::parse();

    let mut pileup = rpileup::PileupIterator::new(&args.input, None, None)?;
    let mut ret: rpileup::IterResult;

    loop {
        println! {"initializing to ref..."}
        ret = pileup.init_to_ref()?;

        match ret {
            rpileup::IterResult::NoData => break,
            _ => loop {
                match pileup.next()? {
                    rpileup::IterResult::ReferenceEnd => break,
                    rpileup::IterResult::NoData => panic!(),
                    _ => (),
                }
            },
        }
    }

    Ok(())
}
