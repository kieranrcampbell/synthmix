
# Synthmix: tools for creating synthetic bulk RNA-sequencing expression
# estimates from single-cell data, using Python, Snakemake + kallisto

# kieran.campbell@sjc.ox.ac.uk
# Copyright (c) 2015

import glob, os
import numpy as np
from synthmix.cm import CreateMix
import json
import sys

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


# Set default config parameters --------------
if not "sc_output_dir" in config.keys():
	config["sc_output_dir"] = os.path.join(config["base_dir"], "scoutput")
if not "bulk_build" in config.keys():
	config["bulk_build"] = os.path.join(config["base_dir"], "bulkbuild")
if not "bulk_output_dir" in config.keys():
	config["bulk_output_dir"] = os.path.join(config["base_dir"], "boutput")
if not "uniform_over_celltypes" in config.keys():
	config["uniform_over_celltypes"] = True
if not "bulk_proportion" in config.keys():
	config["bulk_proportion"] = 0.5
if not "seed" in config.keys():
	config["seed"] = 123
if not "summary_sheet" in config.keys():
	config["summary_sheet"] = "summary_sheet.json"
if not "kallisto" in config.keys():
	config["kallisto"] = "kallisto"
if not "replications" in config.keys():
	config["replications"] = 1

summary_sheet = os.path.join(config["base_dir"], "summary_sheet.json")

# Bulding bulks from single cell ---------------

kallisto_output_files = ['abundance.h5', 'abundance.txt', 'run_info.json']

# sanitise depth input
njobs = len(config["mix_ratio"])
if type(config["depth"]) is list:
	if len(config["depth"]) == 1:
		config["depth"] = [config["depth"][0] for i in range(njobs)]
	else:
		if len(config["depth"] != njobs):
			print("Length of config[\"depth\"] must be 1 or njobs")
			sys.exit(1)
elif type(config["depth"]) is int or type(config["depth"]) is float:
	config["depth"] == [config["depth"] for i in range(njobs)]
else:
	print("config[\"depth\"] must either be a list or int/float")
	sys.exit(1)

replications = config["replications"]

bulk_names = ["" for i in range(replications * njobs)]
bulk_fastq_files = ["" for i in range(replications * njobs)] ## complete list of all files for bulk

if "prefix" in config.keys():
	bulk_names = [config["prefix"] + "_" for i in range(replications * njobs)]

for i in range(njobs):
	for j in range(replications):
		ind = replications * i + j
		bulk_names[ind] += str(config["depth"][i]) + "_"
		bulk_names[ind] += '_'.join([str(m) for m in config["mix_ratio"][i]]).replace('.','')
		bulk_names[ind] += "_id" + str(ind)
		b = os.path.join(config["bulk_build"], "bulk" + bulk_names[ind] + ".fastq.gz")
		bulk_fastq_files[ind] = [b.replace(".fastq", "_" + str(i + 1) + ".fastq") for i in range(2)]

# complete list of all single-cell fastqs - output of build bulk and input to sc quant

scdirs = [os.path.join(config["base_dir"], sc) for sc in config["scdirs"]] # location of single-cell fastqs, also used in quant
cells_full_path = [f for dir in scdirs for f in glob.glob(dir + '/*_*.fastq.gz')]

# Bulk Quantification --------------- 

bulk_quant_dirs = [os.path.join(config["bulk_output_dir"], bn) for bn in bulk_names] # bulk estimate directory
bulk_quant_files = [os.path.join(bqd, kof) for kof in kallisto_output_files for bqd in bulk_quant_dirs] # bulk estimate files


# SIngle-cell quantification ---------------

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

#Pick which cells are being taken forward for bulk comprehension ----- 
# use config['bulk_proportion'] for cell proportion
cells_per_type = [len(x) for x in cells_fastq_path]
cells_per_type_to_choose = [int(round(config['bulk_proportion'] * x)) for x in cells_per_type]

np.random.seed(config["seed"])
cells_for_bulk_index = [np.random.choice(cells_per_type[i],
	cells_per_type_to_choose[i], replace = False) for i in np.arange(len(cells_per_type))]

cells_for_bulk_fastq_path = [] # what we're going to hand to cm.py
for i in range(len(cells_per_type)):
	cells_for_bulk_fastq_path.append([cells_fastq_path[i][j] for j in cells_for_bulk_index[i]])

cells_for_bulk_dict = {config["scdirs"][i]: cells_for_bulk_fastq_path[i] for i in range(len(cells_per_type))}
f = open(summary_sheet, "w")
f.write(json.dumps(cells_for_bulk_dict))
f.close()

# flatten out structure
shell_output_file_list = [o for dir in odirs_cell_level for o in dir]
fastq_strand_1 = [c[0] for type in cells_fastq_path for c in type]
fastq_strand_2 = [c[1] for type in cells_fastq_path for c in type]

#----- Snakemake stuff starts here

INDEX = config["kallisto_index"]
kallisto = config["kallisto"]

# Depth and mixing ratios for each replicate -----

depth = [d for d in config["depth"] for j in range(replications)]
mix_ratio = [m for m in config["mix_ratio"] for j in range(replications)]


# Now for actual snakemake.....
rule buildkallistoindex:
	input:
		config["transcript_index"]
	output:
		INDEX
	shell:
		"{kallisto} index -i {output} {input}"


rule buildbulk:
	input:
		cells_full_path

	output:
		bulk_fastq_files

	run:
		cm = CreateMix(cells_for_bulk_fastq_path, mix_ratio, depth, 
			bulk_fastq_files, paired_end = True, 
			uniform_over_celltypes = config["uniform_over_celltypes"])
		cm.cellIO()

rule quantifysc:
	input:
		cells_full_path
	output:
		sc_all_output_files
	run:
		for output, s1, s2 in zip(shell_output_file_list, fastq_strand_1, fastq_strand_2):
			shell("{kallisto} quant -i {INDEX} -o {output} {s1} {s2}")


rule quantifybulk:
	input:
		bulk_fastq_files
	output:
		bulk_quant_files
	run:
		for job in range(njobs * replications):
			bqd = bulk_quant_dirs[job]
			bff = bulk_fastq_files[job]
			shell("{kallisto} quant -i {INDEX} -o {bqd} {bff[0]} {bff[1]}")

rule all:
	input:
		sc_all_output_files,
		bulk_quant_files
		