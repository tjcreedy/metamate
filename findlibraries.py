#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for matching library-wise sequences with master sequence list and summing reads"""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import Counter
import os

# Global variables

# Class definitions

# Function definitions

# individual libraries read in as sets of sequences - names are not relevant

def match_libraries_to_sequences(master, library_paths):
	"Work through libraries counting incidences of each master sequence"
	
	# Convert master into a simple string dictionary
	master_dict = {str(record.seq).upper() : name for name, record in master.items()}
	
	# TODO: ensure number of items in master_dict matches number of items in master
	
	# Set up empty dictionary for librarywise results
	counts_by_library = dict()
	
	# Set up empty list for all zotu incidences
	all_sequences = []
	
	# Loop through libraries
	for libpath in library_paths:
		
		# Extract library name
		libname = os.path.splitext(os.path.basename(libpath))[0]
		
		# Get list of all zotus in library
		with open(libpath) as fh:
			library_list = [master_dict[seq.upper()] for head, seq in SimpleFastaParser(fh) if seq.upper() in master_dict.keys()]
		
		# Add to librarywise results
		counts_by_library[libname] = Counter(library_list)
		
		# Add to master zotu incidences
		all_sequences.extend(library_list)
	
	# Set up empty dictionary for totals results
	counts_total = { "total" : Counter(all_sequences)}
	
	# Return final dictionary
	return counts_by_library, counts_total




