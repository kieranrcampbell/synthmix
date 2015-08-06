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

	def __init__(self, cell_dirs, mixing_coefficients, 
		read_depth, output_file, paired_end = True,
		uniform_over_celltypes = True):
		""" 
		Construct a new CreateMix

		@param cell_dirs: List of directories holding the fastq files, where
		each directory corresponds to one 'type' of file
		@param mixing_coefficients: A list or tuple where each entry is the 
		proportion of reads corresponding to that type
		@read_depth The read depth of the output file
		@paired_end Is the input & output data paired end?
		@uniform_over_celltypes Logical. If True, the number of reads per cell type
		is taken uniformly across all fastq files for that cell type, so if a given cell has
		more reads then it will be sampled proportional to the number of reads. If 
		False, then a set number of reads is sampled from each cell by taking
		read_depth * mixing_coefficient / n_cells_of_type per cell.
		"""

		opts = { 'cell_dirs': cell_dirs, 
		'mixing_coefficients': mixing_coefficients, 
		'read_depth': read_depth, 
		'output_file': output_file,
		'paired_end': paired_end, 
		'uniform_over_celltypes': uniform_over_celltypes }
		
		self._check_inputs(cell_dirs, mixing_coefficients)

		self.cell_dirs = cell_dirs

		## WARNING: mixing_coefficients explicitly in order of cell_dirs, so
		## if iterating over, make sure keys are taken from cell_dirs and NOT
		## some_dict.keys()
		self.mixing_coefficients = mixing_coefficients
		self.read_depth = read_depth
		self.output_file = output_file
		self.paired_end = paired_end

		""" Dictionary holding list of fastq files for each cell type """
		self.cell_fastq_dict = {}

		""" Dictionary holding list of cell *names* for each cell type """
		self.cell_dict = {}

		""" Dictionary holding list of cell sizes for each type of cell """
		self.cell_sizes = {}

		""" Dictionary holding number of reads coming from each cell """
		self.reads_per_cell = {}


		self.cells_per_type = {}
		

		if self.paired_end:
			self._find_paired_ends()
			self.output_file = [self.output_file.replace(".fastq.gz","") + x + ".fastq.gz" for x in ('_1','_2')]
		else:
			self._trim_dict()
			self.output_file = [self.output_file]

		print('Finding number of cells')
		self._get_ncells()
		print('Finding paired ends')
		self._find_paired_ends()
		print('Calculating reads per cell')
		self._find_reads_per_cell(uniform_over_celltypes)
		print('Init done')

	def _check_inputs(self, cell_directories, mixing_coefficients):

		for d in cell_directories:
			assert os.path.exists(d), "Directory %s does not exit" %d

		assert sum(mixing_coefficients) == 1, "Mixing coefficients must sum to 1"

		assert len(mixing_coefficients) == len(cell_directories), "Must have same size list of selected files as directories: one cell directory per file type"

		return		


	def _get_ncells(self):
		""" Crawls each directory counting the number of fastq files to work 
		out how many reads per cell we should take (at random) """

		for directory in self.cell_dirs:
			self.cell_fastq_dict[directory] = [f for f in os.listdir(directory) if f.endswith(".fastq.gz")]
			self.cells_per_type[directory] = len(self.cell_fastq_dict[directory])

		if self.paired_end:
			assert np.any(np.array(list(self.cells_per_type.values())) % 2 == 0), "Paired end reads must have even number of fastq.gz files"
			for(dir, n_cells) in self.cells_per_type.items():
				self.cells_per_type[dir] = n_cells / 2			


	def  _find_paired_ends(self):
		""" Iterates over self.cell_fastq_dict and sorts files into pairs based on _1 and _2 """

		for celltype, fastq_list in self.cell_fastq_dict.items():
			forward_strand = [x for x in fastq_list if "_1" in x]
			reverse_strand = [x for x in fastq_list if "_2" in x]

			forward_strand_names = sorted([x.replace("_1.fastq.gz","") for x in forward_strand])
			reverse_strand_names = sorted([x.replace("_2.fastq.gz","") for x in reverse_strand])

			assert(forward_strand_names == reverse_strand_names), "Paired end reads must have matching names"
			self.cell_dict[celltype] = forward_strand_names			


	def __countlines(self, path):
		print("Counting lines of " + path)
		nLines = 0
		with gzip.open(path) as gh:
			for line in gh:
				nLines = nLines + 1
		return nLines

	def _find_reads_per_cell(self, uniform_over_celltypes):
		D = self.read_depth

		"""
		For mixing ratios g1, g2, ..., gn, we want g1D, g2D etc
		to come from each type. However, this should be drawn uniformly from
		each cell type, and since each cell has a slightly different read depth, not necessarily
		each cell 
		"""

		reads_per_cell_type = D * np.array(self.mixing_coefficients)

		if uniform_over_celltypes:
			"""
			then reads for cell i given R_j reads for cell type j is
			s_{ij} * R_j / (\sum_j s_{ij}) 
			where s_{ij} is number of reads in cell {ij}, and
			R_j = g_j D
			"""

			# need to count lines in fastq files
			for (dir, cells) in self.cell_dict.items():
				sizes = []
				for cell in cells:
					if self.paired_end:
						cell = cell + "_1.fastq.gz"
					else:
						cell = cell + ".fastq.gz"

					sizes.append(self.__countlines(os.path.join(dir, cell)))
				self.cell_sizes[dir] = sizes

			
			Rj = D * np.array(self.mixing_coefficients)
			for (j, directory) in enumerate(self.cell_dirs):
				s = np.array(self.cell_sizes[directory])
				sprop = s / float(s.sum())
				self.reads_per_cell[directory] = np.round(Rj[j] * sprop)


		else:
			""" 
			reads for cell i given R_j reads for cell type j is
			g_j D / n_j
			for n_j cells of type j
			"""
			for (i, directory) in enumerate(self.cell_dirs):
				self.reads_per_cell[directory] = np.ones(self.cells_per_type[directory]) * \
				reads_per_cell_type[i]  / \
				float(self.cells_per_type[directory])



	def cellIO(self):

		outfilestreams = [gzip.open(f,'wb') for f in self.output_file]

		for cell_type in self.cell_dirs:
			""" iterate over each type of cell """
			## reads_per_cell = self.assigned_reads_per_cell[i] # select reads_per_cell from 
			for (i, cell) in enumerate(self.cell_dict[cell_type]):
				reads_per_cell = self.reads_per_cell[cell_type][i]
				self._write_one_cell(cell_type, cell, reads_per_cell, outfilestreams)


		[f.close() for f in outfilestreams]

	def _write_one_cell(self, directory, cell_file, reads_per_cell, outfilestream):
		""" select reads_per_cell reads at random from cell file in directory, and output
		to outfilestream """
		if self.paired_end:
			#print "Reading %s" % cell_file
			fs = [gzip.open(os.path.join(directory, f),'rb') for f in [cell_file + x + ".fastq.gz" for x in ("_1","_2")]]
			
			lines = [f.readlines() for f in fs]
			n_lines = [len(l) for l in lines]

			assert n_lines[0] == n_lines[1], "Paired end reads must have equal number of lines in _1 and _2 files"
			n_lines = n_lines[0]
			assert n_lines % 4 == 0, "Number of lines in %s  and %s is not a multiple of 4" % cell_file
			

			""" now select reads_per_cell from n_lines """

			chosen_lines = np.random.choice(int(n_lines / 4),  int(reads_per_cell)) * 4

			for l in chosen_lines:
				[outfilestream[i].writelines(lines[i][l:(l+4)]) for i in (0,1)]
			lines = None
			[f.close() for f in fs]

		else:
			print("Non-paired-end currently not supported")
			sys.exit(1)


	def print_all(self):
		print(self.cell_dict)
		print(self.cell_fastq_dict)
		print(self.cell_sizes)
		print(self.reads_per_cell)
		print(self.cells_per_type)




