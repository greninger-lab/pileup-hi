# pileup-hi

<p align="center"> <img src="imgs/pileuphi.png" alt="Peter the parallel
pileup jellyfish" width="65%"/> </p>

<div align="center"> Extensible, high-throughput pileup generation from
SAM/BAM/CRAM <br>
<a href=https://academic.oup.com/bioinformatics/article/42/10/btag721/8842597?login=true#579364384>Manuscript</a>
</div>   

## Table of contents
* [what is pileup\-hi?](#what-is-pileup-hi)
* [what is pileup\-hi <em>not</em>?](#what-is-pileup-hi-not)
* [Installation and demo instructions](#installation-and-demo-instructions)
    * [Requirements](#requirements)
    * [Option 1: install with cargo (recommended)](#option-1-install-with-cargo-recommended)
    * [Option 2: compile from source](#option-2-compile-from-source)
    * [Test command](#test-command)
* [When should I use it over <em>samtools mpileup</em>?](#when-should-i-use-it-over-samtools-mpileup)
* [How can I use it?](#how-can-i-use-it)
    * [Options](#options)
* [Testing](#testing)
* [The new histo format](#the-new-histo-format)

## what is pileup-hi?
pileup-hi is a program for generating pileups from the SAM, BAM, and CRAM formats quickly. pileup-hi implements a pileup-engine with:   
- built-in options for sub-reference multi-threading
- support for custom output formats 
- zero-copy Rust API for creating and iterating over pileup data streams
- binary-identical output to samtools mpileup when using the `mpileup` output format with the binary described in this repo (see [testing](#testing))

this repository contains:
- the pileup-hi_bin binary crate, which behaves similarly to _samtools mpileup_.
- the pileup-hi_lib crate, which implements the underlying pileup engine and API.

## what is pileup-hi _not_?
pileup-hi is not an end-to-end variant calling platform: like _samtools mpileup_, it is a low-level SAM/BAM parser that is meant to be built upon with tools that identify variants, assemble genomes, etc. (e.g. iVar). pileup-hi was developed to retrieve raw nucleotide alignment data from datasets with size/depth that is computationally challenging for other software. 

If an end-to-end variant calling pipeline that handles statistics and math is what you need, check out GATK's HaplotypeCaller and Mutect2 utilities, or iVar if you deal with viruses.

## Installation and demo instructions
#### Requirements
- an installation of htslib 1.2.2 or higher on your system. If you have samtools 1.2.2 or higher on your system, you should be good. Run `samtools --version`, which prints the version of htslib, to check.

If you need to install or upgrade htslib on MacOS, see about installing it via [homebrew](https://brew.sh/) for easy upgrading in the future. For Linux, please use your distribution's package manager.

####  Option 1: install with cargo \(recommended\)
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

#### Test command

You can download and test an input file from this repository's `test_data` directory using the code below. 
```
wget https://raw.githubusercontent.com/epiliper/pileup-hi/main/test_data/c1%23pad1.bam
wget https://raw.githubusercontent.com/epiliper/pileup-hi/main/test_data/c1.fa

pileuphi plp -t 1 -ABQ0 c1#pad1.bam -f c1.fa

# if you have samtools mpileup installed, feel free to compare
samtools mpileup -ABQ0 c1#pad1.bam -f c1.fa

# also, check out the 'histo' mode
pileuphi histo -t 1 -ABQ0 c1#pad1.bam -f c1.fa
```

You should see something like this:  

**pileup-hi plp**
```
c1	1	A	9	^!.^!.^!.^!.^!.^!.^!.^!.^!.	~~~~~~~~~
c1	2	A	9	.........-3CCG	~~~~~~~~~
c1	3	C	9	........*	~~~~~~~~~
c1	4	C	9	........-1G*	~~~~~~~~~
c1	5	G	9	....+6GTTAAC.+6*TTAA*.+6GTT***.+6***AAC*+6**TA**-1C*+6GTTAAC-3CGG	~~~~~~~~~
c1	6	C	9	.......**	~~~~~~~~~
c1	7	G	9	........*	~~~~~~~~~
c1	8	G	9	........*	~~~~~~~~~
c1	9	T	9	.........	~~~~~~~~~
c1	10	T	9	.$.$.$.$.$.$.$.$.$	~~~~~~~~~
```

**samtools mpileup**
```
c1	1	A	9	^!.^!.^!.^!.^!.^!.^!.^!.^!.	~~~~~~~~~
c1	2	A	9	.........-3CCG	~~~~~~~~~
c1	3	C	9	........*	~~~~~~~~~
c1	4	C	9	........-1G*	~~~~~~~~~
c1	5	G	9	....+6GTTAAC.+6*TTAA*.+6GTT***.+6***AAC*+6**TA**-1C*+6GTTAAC-3CGG	   ~~~~~~~~~
c1	6	C	9	.......**	~~~~~~~~~
c1	7	G	9	........*	~~~~~~~~~
c1	8	G	9	........*	~~~~~~~~~
c1	9	T	9	.........	~~~~~~~~~
c1	10	T	9	.$.$.$.$.$.$.$.$.$	~~~~~~~~~
```

**pileup-hi histo**
```
c1	1	9	9	0	0	0	0	0	0	[]	     []
c1	2	9	9	0	0	0	0	0	0	[]	     [1CCG]
c1	3	9	0	0	8	0	0	1	0	[]	     []
c1	4	9	0	0	8	0	0	1	0	[]	     [1G]
c1	5	9	0	7	0	0	0	2	0	[2GTTAAC 1*TTAA* 1GTT*** 1***AAC 1**TA**]	[1C 1CGG]
c1	6	9	0	0	7	0	0	2	0	[]	     []
c1	7	9	0	8	0	0	0	1	0	[]	     []
c1	8	9	0	8	0	0	0	1	0	[]	     []
c1	9	9	0	0	0	9	0	0	0	[]	     []
c1	10	9	0	0	0	9	0	0	0	[]	     []
```


You can also clone this repository and run the entire regression test suite. See [the testing section](#testing).

## When should I use it over _samtools mpileup_?
This question is answered in detail in the manuscript associated with this software, which you can find [here](https://academic.oup.com/bioinformatics/article/42/10/btag721/8842597#579364384).
In short:  

1. You want to process alignments to very long (e.g. eukaryotic) genomes and want to save time.
    - note that parallel execution can work at the sub-reference level, so you don't need a multi-chromosome alignment to use this.
2. You have very high-depth datasets and you wish to save space with the abbreviated histo format, as well as save time.

## How can I use it?
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
- `--rf`: include reads with any of these flags. Accepts a decimal, hexadecimal, or octal SAM flag mask, or comma-separated flag names (for example, `16`, `0x10`, or `READ2,REVERSE`).
- `--ff`: exclude reads with any of these flags. Accepts the same formats as `--rf`.
- `q/--min-MQ`: minimum mapping quality for a  read's bases to be counted. Default = 0.
- `-Q/--min-BQ`: minimum phred score for a base to be counted. Default = 13.
- `-B/--no-BAQ`: Disable BAQ realignment (it's disabled if you don't provide a reference FASTA with `-f`).
- `E/--redo-BAQ`: recalculate BAQ realignment. Requires a reference FASTA.

## Testing
pileup-hi uses the same regression test suite as samtools mpileup, but truncated to the flags pileup-hi supports. You can run tests by cloning this repository and: 

```bash
cd test_data
sh regression.sh mpileup.reg
```

Expect regression tests to be added as more functionality is added.

## The new `histo` format
`histo` was developed as a condensed version of the pileup output format that remains a constant size regardless of alignment depth: it does so by providing a *count* of nucleotides rather than a *list*, and excludes qualities altogether. 

An example of the differences between `plp` and `histo` can be seen in the [Installation and demo instructions section](#test-command).

Below is an explanation of the format:

- column 1: reference name
- column 2: reference coordinate
- column 3: depth (**not including gaps**)
- column 4: #A's
- column 5: #G's
- column 6: #C's
- column 7: #T's
- column 8: #N's
- column 9: #Gaps
- column 10: #Refskips
- Column 11: space-delimited unique insertions\*
- Column 12: space-delimited unique deletions\*

\*Insertions and deletions are represented by their *count* followed by their *sequence*. For example:
```
[4AT 2A] in deletion column (13): 
    - 4 reads w/ deletions of ref bases AT
    - 2 reads w/ deletions of ref base A.
```
