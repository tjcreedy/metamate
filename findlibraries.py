#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for matching library-wise sequences with master sequence list and summing reads"""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter
import os, sys

# Global variables

# Class definitions

# Function definitions

# detect input format
def detect_format(path):
	with open(path) as f:
		fline = f.readline().strip()
	
	if(fline[0] == "@"):
		return("fastq")
	elif(fline[0] == ">"):
		return("fasta")
	else:
		return("unknown")


# individual libraries read in as sets of sequences - names are not relevant
def match_libraries_to_sequences(master, library_paths):
	# master, library_paths = [rawasvs, args.libraries]
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
		#libpath = library_paths[0]
		
		# Extract library name
		libname = os.path.splitext(os.path.basename(libpath))[0]
		
		# Get sequences
		seqformat = detect_format(libpath)
		sequences = None
		with open(libpath) as fh:
			if(seqformat == "fasta"):
				sequences = {head:seq for head,seq in SimpleFastaParser(fh)}
			elif(seqformat == "fastq"):
				sequences = {head:seq for head, seq, qual in FastqGeneralIterator(fh)}
			else:
				sys.stderr.write("Error: can't detect format of " + libpath + "\n")
				exit()
		
		# Get list of all zotus in library
		library_list = [master_dict[seq.upper()] for head, seq in sequences.items() if seq.upper() in master_dict.keys()]
		
		# Add to librarywise results
		counts_by_library[libname] = Counter(library_list)
		
		# Add to master zotu incidences
		all_sequences.extend(library_list)
	
	# Set up dictionary for totals results
	counts_total = { "total" : Counter(all_sequences)}
	
	# Return final dictionary
	return counts_by_library, counts_total
