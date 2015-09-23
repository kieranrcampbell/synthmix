# Synthmix
Create synthetic bulk RNA-seq data from a known mixture of single-cell RNA-seq data, and quantify using [Kallisto](http://pachterlab.github.io/kallisto/).

#### Requirements:
* [Kallisto](http://pachterlab.github.io/kallisto/)
* [Snakemake](https://bitbucket.org/johanneskoester/snakemake/wiki/Home)

## Description

Synthmix takes single-cell fastq files of different cell types and mixes them in a known proportion to a certain read depth in order to mimic bulk RNA-seq data that is a mix of cell types. All fastq files are then quantified to TPM using Kallisto.

Synthmix uses Snakemake for workflow control in order to avoid re-building and re-quantifying files that already exist, which requires Python 3. Note: Synthmix currently only supports paired-end reads in `.fastq.gz` format.

Synthmix is designed so that only a proportion of cells are used to create the bulk data while the rest are quantified at the single-cell level. If all cells were used and quantified then typically predictive algorithms would under-report the error rate, since you can effectively 'see inside' the bulk sample.

## Quick start

To begin, make sure Kallisto is in your path (on linux this can be acheived by `export PATH=/path/to/kallisto/:$PATH`).

Next, you need to tell Synthmix where your files are and what parameters you'd like used to construct the bulk files. This is done via a config file in JSON format:

```json
{
	"base_dir": "/path/to/project",
	"scdirs": ["epithelial","stromal"],
	"mix_ratio": [[0.1, 0.9], [0.2, 0.8]],
	"depth": [1000, 1000],
	"transcript_index": "kallisto_files/transcripts.fasta.gz",
	"kallisto_index": "kallisto_files/transcripts.idx",
	"bulk_proportion": 0.5,
	"seed": 123
}
```

This requires the following fields:
* `base_dir`: Parent directory of where the single-cell fastq files are along with where we'll put all subsequent files
* `scdirs`: Child directories of `base_dir`. Each of these must contain paired end `.fastq.gz` files. Each directory corresponds to a particular cell type.
* `mix_ratio`: A numeric list containing the relative abundances of each cell type (where the order matches that of `scdirs`). Note that this must sum to 1 otherwise an error is thrown.
* `depth`: Number of reads in the output bulk file. Should be a list of number of jobs.
* `transcript_index`: The transcriptome index required by Kallisto
* `kallisto_index`: The pseudoalignment index built by Kallisto (this is created as the first step of Synthmix)
* `bulk_proportion`: The proportion of cells to be mixed for bulk fastq generation. Defaults to 0.5. *Optional*
* `seed`: The random seed to be used to sample cells for bulk/single split. This should be set for your analysis to be reproducible. Defaults to 123. *Optional*
* `summary_sheet`: Json saying which files were witheld and which were taken forward for bulk analysis. *Optional*
* `sc_output_dir`: Where single-cell quantification by Kallisto is stored. Defaults to `os.path.join(base_dir, "scoutput")`. *Optional*
* `bulk_build`: Directory for mixed fastq files. Defaults to `os.path.join(base_dir, "bulkbuild")`. *Optional*
* `bulk_output_dir`: Where bulk quantification by Kallisto is stored. Defaults to `os.path.join(base_dir, "boutput")`. *Optional*
* `uniform_over_celltypes`: Logical. If True, the number of reads per cell type is taken uniformly across all fastq files for that cell type, so if a given cell has more reads then it will be sampled proportional to the number of reads. If False, then a set number of reads is sampled from each cell by taking read_depth * mixing_coefficient / n_cells_of_type per cell. *Optional*

Then call 
```bash
snakemake --configfile synthmix_conf.json
```
Note that any option can also be passed to synthmix through snakemake's built in config rules on the command line, e.g.
```bash
snakemake --configfile synthmix_conf.json --config depth=100
```
By default bulk fastqs are built in a directory called `bulkbuild`, while bulk and single- kallisto quantification output is put in the directories `boutput` and `scoutput` respectively. These can
be changed using the config options `bulk_build`, `sc_output_dir` and `bulk_output_dir`.

Bulk fastqs are built in the format a_b1_b2_s.fastq.gz, where a is the read depth, b1 and b2 are the mixing ratios (can go on to any number) and s is the strand (1 or 2), e.g. 1000_02_08_1.fastq.gz for a bulk sample at a depth of 1000 with mixing proportions of 0.2 and 0.8 respectively.

## Advanced

### Options for Synthmix output directories

Synthmix treats bulk output of identical depth and mixing proportions as the same to avoid re-builds. However, you may wish to produce replicates of the same parameters. To acheive this, specify "prefix" in the configuration and add a unique identifier to force a rebuild under a different name, e.g.

```bash
snakemake --configfile synthmix_conf.json --config prefix=build1
snakemake --configfile synthmix_conf.json --config prefix=build2
```
will create two bulk outputs with the same parameters but with build1 and build2 prepended to differentiate between them (the resultant files will be something like `build1_1000_02_08_1.fastq.gz` for read depth 1000, mixing proportions 0.2 and 0.8 and the first strand).

### Options for bulk data synthesis

Synthmix includes the option `uniform_over_celltypes`, which controls how reads are sampled from single-cells. If `uniform_over_celltypes` is `True` (default) then reads are sampled evenly _for a given cell type_ regardless of fastq size. In other words, if we have two cell types with mixing ratios _a_ and _b_ to a desired read depth _D_ then _aD_ reads will be sampled uniformly from the fastqs corresponding to the first cell type and similarly _bD_ from the second.

However, in general a cell of larger size will contribute more mRNA to a bulk expression estimate than one of smaller size. Therefore, we may want to sample reads from larger cells more often than from smaller cells. If `uniform_over_celltypes` is set to `False` then _aD_ and _bD_ reads are still sampled from each cell type, but rather than uniformly in all cells of that cell type they are sampled in proportion to the size of the fastq file in question.