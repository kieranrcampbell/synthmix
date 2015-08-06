
# Synthmix: tools for creating synthetic bulk RNA-sequencing expression
# estimates from single-cell data, using Python, Snakemake + kallisto

# kieran.campbell@sjc.ox.ac.uk
# Copyright (c) 2015

import glob, os
from synthmix.cm import CreateMix

## from http://www.peterbe.com/plog/uniqifiers-benchmark
def unique(seq, idfun=None): 
   # order preserving
   if idfun is None:
       def idfun(x): return x
   seen = {}
   result = []
   for item in seq:
       marker = idfun(item)
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(item)
   return result

## User required fields
KALLISTO_DIR = "/net/isi-scratch/kieran/tools/kallisto_linux-v0.42.1/"
KALLISTO = os.path.join(KALLISTO_DIR, "kallisto")

BASE_DIR = "/net/isi-project/CW010_CAMPBELL_SCNGEN/data/admix/monocle/"

INDEX = "kallisto_files/transcripts.idx"

SCDIRS = ["type2","type3"]
SC_OUTPUT_DIR = os.path.join(BASE_DIR, "scoutput")

BULK_BUILD = os.path.join(BASE_DIR, "bulkbuild")
MIX_RATIO = [.2, .8]

BULK_OUTPUT_DIR = os.path.join(BASE_DIR, "boutput")

DEPTH = 1000

#----- Bulding bulks from single cell

kallisto_output_files = ['abundance.h5', 'abundance.txt', 'run_info.json']

bulk_name = '_'.join([str(m) for m in MIX_RATIO]).replace('.','')
bulk_file_for_createmix = os.path.join(BULK_BUILD, "bulk" + bulk_name + ".fastq.gz")
bulk_files_for_output = [bulk_file_for_createmix.replace(".fastq", "_" + str(i + 1) + ".fastq") for i in range(2)]

#----- Bulk Quantification

bulk_quant_dir = os.path.join(BULK_OUTPUT_DIR, strmixrat) # bulk estimate directory
bulk_quant_files = [os.path.join(bulk_quant_dir, kof) for kof in kallisto_output_files] # bulk estimate files

scdirs = [os.path.join(BASE_DIR, sc) for sc in SCDIRS] # location of single-cell fastqs, also used in quant

# complete list of all single-cell fastqs - output of build bulk and input to sc quant
cells_full_path = [f for DIR in scdirs for f in glob.glob(DIR + '/*_*.fastq.gz')]

#----- SIngle-cell quantification

odirs = [os.path.join(SC_OUTPUT_DIR, sc) for sc in SCDIRS] # output directories for sc quant

# list of lists - inner corresponds to filenames belonging to a given cell type (outer)
cells = [unique([os.path.basename(f).split('_')[0] for f in glob.glob(d + "/*_*.fastq.gz")]) for d in scdirs]

# list^3 of all fastqs, of the form
# - cell type
# -- cell
# --- fastq_1 fastq_2
cells_fastq_path = [ [ [os.path.join(scdirs[i], cells[i][j] + "_" + str(k + 1) + ".fastq.gz") for k in range(2) ] for j in range(len(cells[i]))    ] for i in range(len(scdirs))]

# list of output directories corresponding to each cell types
odirs_cell_level = [[os.path.join(odirs[i], cells[i][j]) for j in range(len(cells[i]))] for i in range(len(odirs))]

all_output_files = [os.path.join(outdir, kof) for celltype in odirs_cell_level \ 
for outdir in celltype for kof in kallisto_output_files] # complete list of all files generated in sc quant

# flatten out structure
shell_output_file_list = [o for dir in odirs_cell_level for o in dir]
fastq_strand_1 = [c[0] for type in cells_fastq_path for c in type]
fastq_strand_2 = [c[1] for type in cells_fastq_path for c in type]

#----- Snakemake stuff starts here

rule buildkallistoindex:
	input:
		"kallisto_files/transcripts.fasta.gz"
	output:
		{INDEX}		
	shell:
		"{KALLISTO} index -i {output} {input}"


rule buildbulk:
	input:
		cells_full_path

	output:
		bulk_files_for_output

	run:
		cm = CreateMix(scdirs, MIX_RATIO, DEPTH, bulk_file_for_createmix)
		cm.cellIO()

rule quantifysc:
	input:
		cells_full_path
	output:
		all_output_files
	run:
		for output, s1, s2 in zip(shell_output_file_list, fastq_strand_1, fastq_strand_1):
			shell("{KALLISTO} quant -i {INDEX} -o {output} {s1} {s2}")


rule quantifybulk:
	input:
		bulk_files_for_output
	output:
		kallisto_output_files
	rule:
		# same as quantifysc