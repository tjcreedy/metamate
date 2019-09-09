#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by the presence of stop codons in translation"""

# Imports
from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
import scipy.stats
import argparse
import os
import sys
import warnings

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for filtering the sequences in a multifasta according to whether their translation contains stop codons. All sequences must have the same reading frame. The reading frame can be supplied if known or is automatically determined. All sequences must use the same translation table, which follows the NCBI numbering convention (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)")

parser.add_argument("input", help = "input file path", metavar = "FASTA")
parser.add_argument("table", help = "the number referring to the translation table to use", metavar = "TABLE")

parser.add_argument("-r","--reading_frame",		help = "coding frame of sequences, if known", type = int, choices = [1,2,3], metavar = "N")

parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-f","--onefile",		help = "rather than outputting two separate files for passing and failing sequences (the default), output sequences in one file, with the specified suffix appended to the header of failing sequences", default = False, metavar = "SUFF")

parser.add_argument("-c","--detection_confidence", 	help = "confidence level (0 < x < 1) for detection of reading frame (default 0.95, usually no need to change)", type = float, default = 0.95, metavar = "N")
parser.add_argument("-m","--detection_minstops",	help = "minimum number of stops to encounter for detection (default 50, may need to decrease for few sequences)", type = int, default = 100, metavar = "N")
# Function definitions

def detect_frame(seq_records, table, pthresh = 0.95, minstops = 100):
	# Check pthresh
	if not 0 < pthresh < 1:
		sys.exit("Error in detect_frame: pthresh must be greater than 0 and less than 1")
	
	# Set defaults
	counts = [0,0,0]
	seq_records_run = []
	
	# Check input type and assign
	if type(seq_records) is dict:
		seq_records_run = list(seq_records.values())
	elif type(seq_records) is list:
		seq_records_run = seq_records
	elif type(seq_records) is str:
		seq_records_run = SeqIO.parse(seq_records, "fasta")
	else:
		sys.exit("Error in detect_frame: type of seq_records is not dict, list or str")
	
	# Do run to count stops
	n = 0
	p = 0
	
	while((sum(counts) < minstops or p < pthresh)):# Check that there are at least minstops counts and the counts are significantly different - if so stop, else continue
		
		# Extract sequence object depending on input type
		seq_record = ""
		if type(seq_records) is str:
			seq_record = next(seq_records_run)
		else:
			seq_record = seq_records_run[n]
			n+=1
		# Compute counts and add to current counts itemwise
		counts = [ x + y for x,y in zip(counts, stopcount(seq_record, table) ) ]	# Count stops in all reading frames and add to current counts
		
		# Compute p-value if enough total counts
		p = (1 - scipy.stats.chisquare(counts)[1])
	
	# Returns frame with minimum counts and p-value
	return counts.index(min(counts))+1

def stopcount(seq_record, table, frame = (1,2,3)):
	
	# Check input types
	run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
	
	# Run counting
	with warnings.catch_warnings():
		warnings.simplefilter('ignore', BiopythonWarning)
		counts = [seq_record.seq[(i-1):].translate(table = table).count("*") for i in run_frame]
	
	# Return string or list depending on length
	if(len(counts) > 1):
		return counts
	else:
		return counts[0]

if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.input))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	
	# Check for bad options
	
	if((args.detection_confidence != 0.95 or args.detection_minstops != 100) and args.reading_frame):
		print("Warning: specifying a detection confidence or detection minstops is useless if the reading frame is known, this will be ignored")
	
	# Find frame
	
	frame = ""
	pvalue = ""
	if(args.reading_frame):
		frame = args.reading_frame
	else:
		frame = detect_frame(args.input, args.table, args.detection_confidence, args.detection_minstops)
		print("Reading frame %i detected " % (frame))
	
	
	# Output depending on options
	passcount = 0
	failcount = 0
	with open(args.input) as infasta:
		
		if not args.onefile :		# The default, two files
			
			with open(os.path.join(args.output_directory, filename + "_transpass.fa"), "w") as passout:
				with open(os.path.join(args.output_directory, filename + "_transfail.fa"), "w") as failout:
					
					for head, seq in SimpleFastaParser(infasta):
						
						if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
							passout.write(">%s\n%s\n" % (head, seq))
							passcount += 1
						else:
							failout.write(">%s\n%s\n" % (head, seq))
							failcount += 1
			
		else:
			
			with open(os.path.join(args.output_directory, filename + "_transfiltered.fa"), "w") as outfasta:
				
				for head, seq in SimpleFastaParser(infasta):
					
					if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
						outfasta.write(">%s\n%s\n" % (head, seq))
						passcount += 1
					else:
						outfasta.write(">%s\n%s\n" % (head + args.onefile, seq))
						failcount += 1
	
	# Print report
	print("%i of %i sequences passed translation filtering" % (passcount, passcount + failcount))

