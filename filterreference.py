#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by presence in reference set of sequences"""

# Imports

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import argparse
import os

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for filtering the sequences in a multifasta according to whether they exactly match sequences in a second, reference multifasta. All sequences are expected to be in the same reading direction, no reverse-complementation is performed. Sequences must exactly match to pass")

parser.add_argument("input", help = "input file path", metavar = "FASTA")
parser.add_argument("reference", help = "reference file path", metavar = "REF")

parser.add_argument("-T", "--threads",			help = "number of threads to use", default = 3, metavar = "N")

parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-f","--onefile",		help = "rather than outputting two separate files for passing and failing sequences (the default), output sequences in one file, with the specified suffix appended to the header of failing sequences", default = False, metavar = "SUFF")

parser.add_argument("-j", "--match_length", help = "the minimum alignment length to consider a BLAST match when comparing zotus against the reference", type = int, metavar = "N", default = 350)
parser.add_argument("-i", "--match_percent", help = "the minimum percent identity to consider a BLAST match when comparing zotus against the reference", type = int, metavar = "N", default = 97)


# Function definitons

def blast_for_presence(queryfile, subjectfile, perc_identity, min_length, threads):
	
	# Set up blast object
	blastn_cline = NcbiblastnCommandline(query = queryfile, subject = subjectfile,
							   evalue = 0.001, perc_identity = perc_identity,
							   num_threads = threads, #TODO set up global threads variable
							   outfmt = 5, out = "refblast.xml")
	
	# Run blast
	stdout, stderr = blastn_cline()
	
	# Get result
	blast_result_fh = open("refblast.xml")
	blast_records = NCBIXML.parse(blast_result_fh)
	
	# Work through results finding reference hits that are suitable
	
	query_pass = set()
	for blast_record in blast_records:
		length_pass = []
		for alignment in blast_record.alignments:
			length_pass.extend([hsp.query_end - hsp.query_start + 1 > min_length for hsp in alignment.hsps])
		if(len(length_pass) > 0 and any(length_pass)):
			query_pass.add(blast_record.query)
	
	return(query_pass)




if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.input))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	
	# Check for bad options
	
	# Do blast search
	
	passed_headers = blast_for_presence(args.input, args.reference, 97, 350, args.threads)
	
	# Output depending on options
	passcount = 0
	failcount = 0
	with open(args.input) as infasta:
		
		if not args.onefile :		# The default, two files
			
			with open(os.path.join(args.output_directory, filename + "_inref.fa"), "w") as passout:
				with open(os.path.join(args.output_directory, filename + "_outref.fa"), "w") as failout:
					
					for head, seq in SimpleFastaParser(infasta):
						
						if head in passed_headers:
							passout.write(">%s\n%s\n" % (head, seq))
							passcount += 1
						else:
							failout.write(">%s\n%s\n" % (head, seq))
							failcount += 1
			
		else:
			
			with open(os.path.join(args.output_directory, filename + "_lengthfiltered.fa"), "w") as outfasta:
				
				for head, seq in SimpleFastaParser(infasta):
					
					if head in passed_headers:
						outfasta.write(">%s\n%s\n" % (head, seq))
						passcount += 1
					else:
						outfasta.write(">%s\n%s\n" % (head + args.onefile, seq))
						failcount += 1
	
	# Print report
	print("%i of %i sequences were in the reference" % (passcount, passcount + failcount))

