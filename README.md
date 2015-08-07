# Synthmix
Create synthetic bulk RNA-seq data from a known mixture of single-cell RNA-seq data, and quantify using [Kallisto](http://pachterlab.github.io/kallisto/).

## Description

Synthmix takes single-cell fastq files of different cell types and mixes them in a known proportion to a certain read depth in order to mimic bulk RNA-seq data that is a mix of cell types. All fastq files are then quantified to TPM using Kallisto.

Synthmix uses Snakemake for workflow control in order to avoid re-building and re-quantifying files that already exist, which requires Python 3. Note: Synthmix currently only supports paired-end reads in `.fastq.gz` format.

## Quick start

To begin, make sure Kallisto is in your path (on linux this can be acheived by `export PATH=/path/to/kallisto/:$PATH`).

Next, you need to tell Synthmix where your files are and what parameters you'd like used to construct the bulk files. This is done via a config file in JSON format:

```json
{
	"base_dir": "/path/to/project",
	"scdirs": ["epithelial","stromal"],
	"mix_ratio": [0.1, 0.9],
	"depth": 1000,
	"transcript_index": "kallisto_files/transcripts.fasta.gz",
	"kallisto_index": "kallisto_files/transcripts.idx"
}
```

This requires the following fields:
* `base_dir`: Parent directory of where the single-cell fastq files are along with where we'll put all subsequent files
* `scdirs`: Child directories of `base_dir`. Each of these must contain paired end `.fastq.gz` files
* `mix_ratio`: A numeric list containing the relative abundances of each cell type (where the order matches that of `scdirs`). Note that this must sum to 1 otherwise an error is thrown.
* `depth`: Number of reads in the output bulk file
* `transcript_index`: The transcriptome index required by Kallisto
* `kallisto_index`: The pseudoalignment index built by Kallisto (this is created as the first step of Synthmix)


Then call 
```bash
snakemake --configfile synthmix_conf.json
```

## Advanced

### Options for Synthmix output directories

### Options for bulk data synthesis