"""

Create synthetic bulk samples from single-cell data


kieran.campbell@sjc.ox.ac.uk

"""

from __future__ import print_function, division
import gzip
import numpy as np
import sys
import os
import logging


class CreateMix:
	""" Take single-cell data and create synthetic bulk data
	"""

	def __init__(self, cell_fastq_list, mixing_coefficients, 
		read_depth, output_files, paired_end = True,
		uniform_over_celltypes = True):
		""" 
		Construct a new CreateMix

		@param cell_fastq_list: A list of depth three: the outer list corresponds
		to different cell *types*, the second list corresponds to different cells,
		and the inner list should always be of length 2 and correspond to the different
		strands (_1 and _2) for synthesis.
		@param mixing_coefficients: A list of lists - each entry corresponds to a 'job' and 
		each 'job' within that will be a list of coefficients
		@read_depth A list with a read depth for each 'job'
		@output_files A list of length(njobs) where each entry is a list of length 2 with
		forward and reverse strand info		
		@paired_end Is the input & output data paired end?
		@uniform_over_celltypes Logical. If True, the number of reads per cell type
		is taken uniformly across all fastq files for that cell type, so if a given cell has
		more reads then it will be sampled proportional to the number of reads. If 
		False, then a set number of reads is sampled from each cell by taking
		read_depth * mixing_coefficient / n_cells_of_type per cell.
		"""

		""" NB order always given by cell_fastq_list """ 

		self.cell_fastq = cell_fastq_list

		## WARNING: mixing_coefficients explicitly in order of cell_dirs, so
		## if iterating over, make sure keys are taken from cell_dirs and NOT
		## some_dict.keys()
		self.mixing_coefficients = mixing_coefficients
		self.read_depth = read_depth
		self.output_files = output_files
		self.paired_end = paired_end
		self.uniform_over_celltypes = uniform_over_celltypes
		self.njobs = len(mixing_coefficients)

		print("[CreateMix] Total jobs: " + str(self.njobs))

		## put some numbers on things
		self.n_cell_types = len(cell_fastq_list)
		self.cells_per_type = [len(x) for x in cell_fastq_list]

		## check everything is in order
		self._check_inputs()


		# """ Dictionary holding list of fastq files for each cell type """
		# self.cell_fastq_dict = {}

		# """ Dictionary holding list of cell *names* for each cell type """
		# self.cell_dict = {}

		""" List holding list of cell sizes for each type of cell """
		self.cell_sizes = []

		""" List holding number of reads coming from each cell *for each job* """
		self.reads_per_cell = []


		# if self.paired_end:
		# 	# self._find_paired_ends()
		# 	for i in range(self.njobs):
		# 		self.output_files[i] = [self.output_files[i].replace(".fastq.gz","") + x + ".fastq.gz" for x in ('_1','_2')]
		# else:
		# 	# self._trim_dict() # what is this supposed to do?
		# 	self.output_file = [self.output_file]

		self.opts = { 'cell_fastq': cell_fastq_list, 
		'mixing_coefficients': mixing_coefficients, 
		'read_depth': read_depth, 
		'output_files': output_files,
		'paired_end': paired_end, 
		'uniform_over_celltypes': uniform_over_celltypes,
		'n_cell_types': self.n_cell_types,
		'cells_per_type': self.cells_per_type,
		'cell_sizes': self.cell_sizes}


		self._find_reads_per_cell()
		

	def _check_inputs(self):

		# for d in cell_directories:
			# assert os.path.exists(d), "Directory %s does not exit" %d

		""" assert number of jobs is consistent """
		assert self.njobs == len(self.output_files), "Must have mixing coefficient and output file for each job "
		assert len(self.read_depth) == self.njobs, "Must have depth for each job" # TODO change this


		for j in range(self.njobs):
			assert sum(self.mixing_coefficients[j]) == 1, "Mixing coefficients must sum to 1"
			assert len(self.mixing_coefficients[j]) == self.n_cell_types, "Must have same size list of selected files as directories: one cell directory per file type"

		self._check_output_consistency()

		return		

	# def _get_cell_name(self, cell_path):
	# 	""" given "/path/to/cellname_x.fastq.gz" return "cellname" """
	# 	s = cell_path.split("/")[-1]
	# 	s = s.replace(".fastq.gz", "")
	# 	s = s.replace("_1", "").replace("_2", "")		
	# 	return s


	# def _get_ncells(self):
	# 	""" Crawls each directory counting the number of fastq files to work 
	# 	out how many reads per cell we should take (at random) """

	# 	for directory in self.cell_dirs:
	# 		self.cell_fastq_dict[directory] = [f for f in os.listdir(directory) if f.endswith(".fastq.gz")]
	# 		self.cells_per_type[directory] = len(self.cell_fastq_dict[directory])

	# 	if self.paired_end:
	# 		assert np.any(np.array(list(self.cells_per_type.values())) % 2 == 0), "Paired end reads must have even number of fastq.gz files"
	# 		for(dir, n_cells) in self.cells_per_type.items():
	# 			self.cells_per_type[dir] = n_cells / 2			


	# def  _find_paired_ends(self):
	# 	""" Iterates over self.cell_fastq_dict and sorts files into pairs based on _1 and _2 """

	# 	for celltype, fastq_list in self.cell_fastq_dict.items():
	# 		forward_strand = [x for x in fastq_list if "_1" in x]
	# 		reverse_strand = [x for x in fastq_list if "_2" in x]

	# 		forward_strand_names = sorted([x.replace("_1.fastq.gz","") for x in forward_strand])
	# 		reverse_strand_names = sorted([x.replace("_2.fastq.gz","") for x in reverse_strand])

	# 		assert(forward_strand_names == reverse_strand_names), "Paired end reads must have matching names"
	# 		self.cell_dict[celltype] = forward_strand_names			


	def __countlines(self, path):
		print("[CreateMix] Counting lines")
		# print("Counting lines of " + path)
		nLines = 0
		with gzip.open(path) as gh:
			for line in gh:
				nLines = nLines + 1
		return nLines

	def _find_reads_per_cell(self):
		

		"""
		For mixing ratios g1, g2, ..., gn, we want g1D, g2D etc
		to come from each type. However, this should be drawn uniformly from
		each cell type, and since each cell has a slightly different read depth, not necessarily
		each cell 
		"""



		if self.uniform_over_celltypes:
			"""
			then reads for cell i given R_j reads for cell type j is
			s_{ij} * R_j / (\sum_j s_{ij}) 
			where s_{ij} is number of reads in cell {ij}, and
			R_j = g_j D
			"""

			# need to count lines in fastq files
			# for (dir, cells) in self.cell_dict.items():
			# 	sizes = []
			# 	for cell in cells:
			# 		if self.paired_end:
			# 			cell = cell + "_1.fastq.gz"
			# 		else:
			# 			cell = cell + ".fastq.gz"

			# 		sizes.append(self.__countlines(os.path.join(dir, cell)))
			# 	self.cell_sizes[dir] = sizes


			# NEW
			for cell_type in self.cell_fastq:
				sizes = []
				for cell in cell_type:
					if self.paired_end: # it's a list!
						fastq_path = cell[0]
					else:
						fastq_path = cell
					sizes.append(self.__countlines(fastq_path))
				self.cell_sizes.append(sizes)

			for job in range(self.njobs):
				D = self.read_depth[job]
				Rj = D * np.array(self.mixing_coefficients[job])
				reads_per_cell_this_job = []
				for j in range(self.n_cell_types):
					s = np.array(self.cell_sizes[j])
					sprop = s / float(s.sum())
					reads_per_cell_this_job.append(np.round(Rj[j] * sprop))
				self.reads_per_cell.append(reads_per_cell_this_job)


		else:
			""" 
			reads for cell i given R_j reads for cell type j is
			g_j D / n_j
			for n_j cells of type j
			"""
			for job in range(self.njobs):
				reads_per_cell_type = D[j] * np.array(self.mixing_coefficients[j])
				reads_per_cell_this_job = []
				for j in range(self.n_cell_types):
					reads_per_cell_this_job.append( np.ones(self.cells_per_type[j]) * \
					reads_per_cell_type[j]  / \
					float(self.cells_per_type[j]) )
				self.reads_per_cell.append(reads_per_cell_this_job)

	def _check_output_consistency(self):
		""" checks _1 & _2.fastq.gz are present in output files correctly """ 
		for job in range(len(self.output_files)):
			assert "_1.fastq.gz" in self.output_files[job][0], "Output missing first strand"
			assert "_2.fastq.gz" in self.output_files[job][1], "Output missing second strand"

	def cellIO(self):
		print("[CreateMix] Writing files")

		outfilestreams = [[gzip.open(f,'wb') for f in job] for job in self.output_files]

		for j in range(len(self.cell_fastq)):  # iterate over each type of cell
			for i in range(len(self.cell_fastq[j])): # iterate over each cell
				""" hand off to read one cell with input files defined in 
				self.cell_fastq[i][j] and a predetermined number of reads nreads
				to one (single) or two (paired-end) output files in outfilestreams """
				self._write_one_cell(j, i, outfilestreams)
				#self._write_one_cell(cell_type, cell, reads_per_cell, outfilestreams)


		for job in outfilestreams:
			for f in job:
				f.close()

		print("[CreateMix] Done")

	def _write_one_cell(self, j, i, outfilestreams):
		""" 
		This function deals with one cell (cell type = j, cell = i) and does IO for each job.
		The logic is we only read each file once saving time.
		"""
		if self.paired_end:
			#print "Reading %s" % cell_file
			fs = [gzip.open(f,'rb') for f in self.cell_fastq[j][i]]
			
			lines = [f.readlines() for f in fs]
			[f.close() for f in fs]
			n_lines = [len(l) for l in lines]

			assert n_lines[0] == n_lines[1], "Paired end reads must have equal number of lines in _1 and _2 files"
			n_lines = n_lines[0]
			assert n_lines % 4 == 0, "Number of lines in %s  and %s is not a multiple of 4" % cell_file
			
			for job in range(self.njobs):
				outfilestream = outfilestreams[job]
				#print(job)
				#print(outfilestream)
				nreads = self.reads_per_cell[job][j][i]
				chosen_lines = np.random.choice(int(n_lines / 4),  int(nreads)) * 4
				# print(chosen_lines)
				for l in chosen_lines:
					[outfilestream[k].writelines(lines[k][l:(l+4)]) for k in (0,1)]
			# lines = None

		else:
			print("Non-paired-end currently not supported")
			sys.exit(1)


	def print_all(self):
		print(self.opts)






