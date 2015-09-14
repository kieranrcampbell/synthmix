
# Synthmix: tools for creating synthetic bulk RNA-sequencing expression
# estimates from single-cell data, using Python, Snakemake + kallisto

# kieran.campbell@sjc.ox.ac.uk
# Copyright (c) 2015

import glob, os
import numpy as np
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

# exec(open("config.py").read()) # read in config


if not "sc_output_dir" in config.keys():
	config["sc_output_dir"] = os.path.join(config["base_dir"], "scoutput")
if not "bulk_build" in config.keys():
	config["bulk_build"] = os.path.join(config["base_dir"], "bulkbuild")
if not "bulk_output_dir" in config.keys():
	config["bulk_output_dir"] = os.path.join(config["base_dir"], "boutput")
if not "uniform_over_celltypes" in config.keys():
	config["uniform_over_celltypes"] = True
if not "sample_mix" in config.keys():
	config["sample_mix"] = [0.5, 0.5]
if not "seed" in config.keys():
	config["seed"] = 123


#--------------- Bulding bulks from single cell

kallisto_output_files = ['abundance.h5', 'abundance.txt', 'run_info.json']

bulk_name = ""
if "prefix" in config.keys():
	bulk_name = config["prefix"]

bulk_name = str(config["depth"]) + "_"
bulk_name += '_'.join([str(m) for m in config["mix_ratio"]]).replace('.','')
bulk_file_for_createmix = os.path.join(config["bulk_build"], "bulk" + bulk_name + ".fastq.gz")
bulk_fastq_files = [bulk_file_for_createmix.replace(".fastq", "_" + str(i + 1) + ".fastq") for i in range(2)]

# complete list of all single-cell fastqs - output of build bulk and input to sc quant

scdirs = [os.path.join(config["base_dir"], sc) for sc in config["scdirs"]] # location of single-cell fastqs, also used in quant
cells_full_path = [f for dir in scdirs for f in glob.glob(dir + '/*_*.fastq.gz')]

#--------------- Bulk Quantification

bulk_quant_dir = os.path.join(config["bulk_output_dir"], bulk_name) # bulk estimate directory
bulk_quant_files = [os.path.join(bulk_quant_dir, kof) for kof in kallisto_output_files] # bulk estimate files


#--------------- SIngle-cell quantification

odirs = [os.path.join(config["sc_output_dir"], sc) for sc in config["scdirs"]] # output directories for sc quant

# list of lists - inner corresponds to filenames belonging to a given cell type (outer)
cells = [unique([os.path.basename(f).split('_')[0] for f in glob.glob(d + "/*_*.fastq.gz")]) for d in scdirs]

# list^3 of all fastqs, of the form
# - cell type
# -- cell
# --- fastq_1 fastq_2
cells_fastq_path = [ [ [os.path.join(scdirs[i], cells[i][j] + "_" + str(k + 1) + ".fastq.gz") for k in range(2) ] for j in range(len(cells[i]))    ] for i in range(len(scdirs))]

# list of output directories corresponding to each cell types
odirs_cell_level = [[os.path.join(odirs[i], cells[i][j]) for j in range(len(cells[i]))] for i in range(len(odirs))]

sc_all_output_files = [os.path.join(outdir, kof) for celltype in odirs_cell_level \ 
for outdir in celltype for kof in kallisto_output_files] # complete list of all files generated in sc quant

# flatten out structure
shell_output_file_list = [o for dir in odirs_cell_level for o in dir]
fastq_strand_1 = [c[0] for type in cells_fastq_path for c in type]
fastq_strand_2 = [c[1] for type in cells_fastq_path for c in type]

#----- Snakemake stuff starts here

INDEX = config["kallisto_index"]

rule buildkallistoindex:
	input:
		config["transcript_index"]
	output:
		INDEX
	shell:
		"kallisto index -i {output} {input}"


rule buildbulk:
	input:
		cells_full_path

	output:
		bulk_fastq_files

	run:
		cm = CreateMix(scdirs, config["mix_ratio"], config["depth"], 
			bulk_file_for_createmix, paired_end = True, 
			uniform_over_celltypes = config["uniform_over_celltypes"])
		cm.cellIO()

rule quantifysc:
	input:
		cells_full_path
	output:
		sc_all_output_files
	run:
		for output, s1, s2 in zip(shell_output_file_list, fastq_strand_1, fastq_strand_2):
			shell("kallisto quant -i {INDEX} -o {output} {s1} {s2}")


rule quantifybulk:
	input:
		bulk_fastq_files
	output:
		bulk_quant_files
	shell:
		"kallisto quant -i {INDEX} -o {bulk_quant_dir} {bulk_fastq_files[0]} {bulk_fastq_files[1]}"

rule all:
	input:
		sc_all_output_files,
		bulk_quant_files
		