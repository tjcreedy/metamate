#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by length variants"""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import argparse
import os

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for filtering the sequences in a multifasta according to their length. Length can be specified by an appropriate combination of the optional arguments: either 1) the expected length, 2) one or both of minimum or maximum length, or 3) a combination of the expected length and some type of variation. Filtering can be strictly limited to sequences varying by a whole number of codons from the expected length")

parser.add_argument("input", help = "input file path", metavar = "FASTA")

parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")

parser.add_argument("-n", "--minimum_length",		help = "remove sequences that are shorter than this value", type = int, default = 0, metavar = "N")
parser.add_argument("-x", "--maximum_length",		help = "remove sequences that are longer than this value", type = int, default = float('Inf'), metavar = "N")
parser.add_argument("-l", "--expected_length",		help = "the expected length of the sequences", type = int, default = 0, metavar = "N")
parser.add_argument("-p", "--percent_variation",	help = "the percentage variation from the expected length within which sequences should be kept", type = float, default = 0, metavar = "N")
parser.add_argument("-b", "--bases_variation",		help = "the number of bases of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")
parser.add_argument("-c", "--codons_variation",		help = "the number of codons of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")

parser.add_argument("-f","--onefile",			help = "rather than outputting two separate files for passing and failing sequences (the default), output sequences in one file, with the specified suffix appended to the header of failing sequences", default = False, metavar = "SUFF")
parser.add_argument("--only_vary_by_codon",		help = "remove sequences that do not vary by a multiple of 3 bases from the expected length", action = "store_true")


# Class definitons

# Function definitions

def resolve_length_spec(length_spec):
	"Resolves specifications for length filtering to just minimum and maximum"
	
	# Set defaults
	min_bp_out = 0
	max_bp_out = float('Inf')
	
	# Compute some useful values for checking specifications
	n_specs = len(length_spec)
	min_max_n = ["min_bp" in length_spec, "max_bp" in length_spec].count(True)
	
	# Check that the specifications match the allowable sets, and if not die
	if not ( n_specs == 0 or										# No specifications at all
		( n_specs == 1 and ( min_max_n == 1 or "len_bp" in length_spec) ) or				# One of minimum, maximum or specified length
		( ( n_specs == 2 and  min_max_n == 2 ) or ( "len_bp" in length_spec and min_max_n == 0 ) )	# Either minimum and maximum, or specified length and a variation around it
	       ) :
		sys.exit("Error: insufficient or unclear length specification")
	
	if n_specs == 0:
		print("Warning: no length specification supplied, all sequences will pass filtering")
	
	# Overwrite min or max values if supplied
	if "min_bp" in length_spec:
		min_bp_out = length_spec["min_bp"]
	if "max_bp" in length_spec:
		max_bp_out = length_spec["max_bp"]
	
	# Set min or max values if length specified
	if "len_bp" in length_spec:
		# Set default
		var_bp = 0
		
		# Overwrite with (calculated) value if variation specified
		if "var_bp" in length_spec:
			var_bp = length_spec["var_bp"]
		elif "var_pc" in length_spec:
			var_bp = round( (length_spec["var_pc"]/100) * length_spec["len_bp"] , 0)
		elif "var_cdn" in length_spec:
			var_bp = length_spec["var_cdn"] * 3
		
		# Overwrite min or max values
		min_bp_out = length_spec["len_bp"] - var_bp
		max_bp_out = length_spec["len_bp"] + var_bp
	
	# Return final values
	return (min_bp_out, max_bp_out)

def check_length(seq_record, min_max_bp, len_bp, var_by_codon = False):
	"Checks if a sequence record (or sequence string) passes or fails the length specifications"
	
	# Check length is ok
	if len_bp == 0 and var_by_codon:
		sys.exit("Error: check_length requires len_bp > 0 if var_by_codon == True")
	
	# Find length
	curr_len = len(seq_record)
	
	# Determine if length falls outside specification
	return ( curr_len >= min_max_bp[0] and curr_len <= min_max_bp[1] ) and ( not var_by_codon or abs(len_bp - curr_len) % 3 == 0 )


if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.input))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	
	# Check for bad options
	
	if args.only_vary_by_codon and not args.expected_length > 0:
		sys.exit("Expected length must be specified if strict codon variation is desired")
	
	# Create dictionary of length specifications
	
	length_spec = {}
	if args.minimum_length > 0:		length_spec["min_bp"] = args.minimum_length
	if args.maximum_length < float('Inf'):	length_spec["max_bp"] = args.maximum_length
	if args.expected_length > 0:		length_spec["len_bp"] = args.expected_length
	if args.bases_variation > 0:		length_spec["var_bp"] = args.bases_variation
	if args.percent_variation > 0:		length_spec["var_pc"] = args.percent_variation
	if args.percent_variation > 0:		length_spec["var_cdn"] = args.codons_variation
	
	# Resolve length specifications
	
	min_max_bp = resolve_length_spec(length_spec)
	
	# Output depending on options
	passcount = 0
	failcount = 0
	with open(args.input) as infasta:
		
		if not args.onefile :		# The default, two files
			
			with open(os.path.join(args.output_directory, filename + "_lengthpass.fa"), "w") as passout:
				with open(os.path.join(args.output_directory, filename + "_lengthfail.fa"), "w") as failout:
					
					for head, seq in SimpleFastaParser(infasta):
						
						if check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) :
							passout.write(">%s\n%s\n" % (head, seq))
							passcount += 1
						else:
							failout.write(">%s\n%s\n" % (head, seq))
							failcount += 1
			
		else:
			
			with open(os.path.join(args.output_directory, filename + "_lengthfiltered.fa"), "w") as outfasta:
				
				for head, seq in SimpleFastaParser(infasta):
					
					if check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) :
						outfasta.write(">%s\n%s\n" % (head, seq))
						passcount += 1
					else:
						outfasta.write(">%s\n%s\n" % (head + args.onefile, seq))
						failcount += 1
	
	# Print report
	print("%i of %i sequences passed length filtering" % (passcount, passcount + failcount))


