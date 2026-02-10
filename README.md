# pileup-hi
## what is pileup-hi, and what is it not? 
pileup-hi is a high-throughput pileup engine for the SAM/BAM file formats. It is multi-threaded and supports the development of custom output formats. This repository contains code to compile the CLI program itself, as well as library code detailing how to construct your own output format (`src/output.rs`).

pileup-hi is not an end-to-end variant calling platform: like _samtools mpileup_, it is a low-level SAM/BAM parser that is meant to be built upon with tools that identify variants, assemble genomes, etc. (e.g. iVar). pileup-hi was developed to retrieve raw nucleotide alignment data from datasets with size/depth that is computationally challenging for other software. 

If an end-to-end variant calling pipeline that handles statistics and math is what you need, check out GATK's HaplotypeCaller and Mutect2 utilities, or iVar if you deal with viruses.

pileup-hi currently has two subcommands that dictate which output format it emits: 

- `plp`: standard _samtools mpileup_ output format.
- `histo`: a list of nucleotide and indel frequencies per coordinate, essentially a condensed form of the pileup format that doesn't grow linearly with alignment depth.

## how can I install it? 
#### Requirements
- an installation of htslib 1.2.2 or higher on your system. If you have samtools 1.2.2 or higher on your system, you should be good. Run `samtools --version`, which prints the version of htslib, to check.

If you need to install or upgrade htslib on MacOS, see about installing it via [homebrew](https://brew.sh/) for easy upgrading in the future. For Linux, please use your distribution's package manager.

####  Option 1. install with cargo
If you don't have cargo installed, see [here](https://rust-lang.org/tools/install/).
```bash
cargo install pileup-hi
```

#### Option 2: compile from source
```bash
git clone https://github.com/epiliper/pileup-hi
cd pileup-hi
make build
```

## when should I use it over _samtools mpileup_?
This question is answered in detail in the manuscript associated with this software, which you can find [here]().
In short:  

1. You want to process alignments to very long (e.g. eukaryotic) genomes and want to save time.
2. You have very high-depth datasets and you wish to save space with the abbreviated histo format, as well as save time.
3. When you want to leverage multiple CPU cores to gain a speedup with longer genomes (even if your BAM only has one reference).

## how can I use it?
pileup-hi attempts to stay consistent with the CLI of samtools mpileup when possible. Below are a list of parameters you can supply. Run `pileuphi <plp|histo> --help` for more information.

usage: `pileuphi <COMMAND> <FILE> [OPTIONS]`

where `COMMAND` is either `histo` or `plp` for the different output formats. Other arguments are shared for both commands.

`FILE` is required and can be either a BAM/SAM file or `stdout`.

### Options
- `-r/--region`: Restrict analysis to a region of a reference. Format of this argument is the same as in samtools mpileup: `REFNAME:START-END`, e.g. `Chr1:400-800` for coordinates 400-800, `Chr1` for all coordinates for Chr1, or `Chr1:400` for all coordinates of Chr1 after and including 400.

- `-o/--output`: file for output. default = `stdout`
- `-a`: emit output for coordinates with no coverage if reference has coverage anywhere.
- `--a`: emit output for all coordinates of all references, regardless of whether or not they have coverage.
- `-t/--threads`: number of threads to use. Requires that the input file be indexed. One thread can process a minimum of 250,000 reference coordinates, though this can be adjusted with `-c/--thread-coords`. Default: number of CPU cores.
- `-c/--thread-coords`: The minimum number of reference coordinates that can be processed by a single worker thread (how small the per-thread chunks can be). 250,000 is an empirically sane value that prevents a lot of repeat work in my experience, and adjusting is not recommended.
- `f/--fasta-ref`: Path to FASTA reference file. 
- `-d/--depth`: Restrict read buffer to only hold a maximum X number of reads per coordinate. Passing 0 is equivalent to specifying infinite depth. This option will lead to output differences between single- and multi-threaded modes.
- `x`: disable R1/R2 overlap correction.
- `A`: count reads with unmapped mates.
- `--rf`: include reads with these flags only.
- `--ff`: exclude reads with these flags.
- `q/--min-MQ`: minimum mapping quality for a  read's bases to be counted. Default = 0.
- `-Q/--min-BQ`: minimum phred score for a base to be counted. Default = 13.
- `-B/--no-BAQ`: Disable BAQ realignment (it's disabled if you don't provide a reference FASTA with `-f`).
- `E/--redo-BAQ`: recalculate BAQ realignment. Requires a reference FASTA.
