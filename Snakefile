
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

KALLISTO_DIR = "/net/isi-scratch/kieran/tools/kallisto_linux-v0.42.1/"
KALLISTO = os.path.join(KALLISTO_DIR, "kallisto")

INDEX = "kallisto_files/transcripts.idx"

BASE_DIR = "/net/isi-project/CW010_CAMPBELL_SCNGEN/data/admix/monocle/"
SCDIRS = ["type2","type3"]
OUTPUT_DIR = BASE_DIR + "output"

scdirs = [os.path.join(BASE_DIR, sc) for sc in SCDIRS]
odirs = [os.path.join(OUTPUT_DIR, sc) for sc in SCDIRS]

# cells = [os.path.basename(f)  for DIR in scdirs for f in glob.glob(DIR + '/*_*.fastq.gz')]
cells_full_path = [f for DIR in scdirs for f in glob.glob(DIR + '/*_*.fastq.gz')]

cells = [unique([os.path.basename(f).split('_')[0] for f in glob.glob(d + "/*_*.fastq.gz")]) for d in scdirs]

cells_fastq_path = [ [ [os.path.join(scdirs[i], cells[i][j] + "_" + str(k + 1) + ".fastq.gz") for k in range(2) ] for j in range(len(cells[i]))    ] for i in range(len(scdirs))]

odirs_cell_level = [[os.path.join(odirs[i], cells[i][j]) for j in range(len(cells[i]))] for i in range(len(odirs))]

kallisto_output_files = ['abundance.h5', 'abundance.txt', 'run_info.json']

all_output_files = [os.path.join(outdir, kof) for celltype in odirs_cell_level for outdir in celltype for kof in kallisto_output_files]

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
		# all files in single cell directories	

	output:
		# bulk files

	shell:
		"python -m synthmix ALL_PARAMS"

rule quantifysc:
	input:
		cells_full_path
	output:
		all_output_files
	run:
		for output, s1, s2 in zip(shell_output_file_list, fastq_strand_1, fastq_strand_1):
			shell("{KALLISTO} quant -i {INDEX} -o {output} {s1} {s2}")


# rule quantifybulk:
# 	input:
# 		# bulk files made by buildbulk
# 	output:
# 		# kallisto files 
# 	rule:
# 		# same as quantifysc