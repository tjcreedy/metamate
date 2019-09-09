#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" """

# Imports

import sys
import argparse
import os

from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

import re

import filterlength
import filterreference
import filtertranslate

import findlibraries
import findclades
import findtaxa

import categorycounting
import assessmentcore


# Global variables

parser = argparse.ArgumentParser(description = "")

	# Input / output paths
parser.add_argument("-Z", "--zotus",			help = "path to a fasta of unique sequences to filter", metavar = "ZOTUs", type = str, default = "test/bee_mocks/all.zotus.glopairmaxit1k_align.fa")
parser.add_argument("-R", "--reference",		help = "path to a fasta of known correct reference sequences", metavar = "REF", type = str, default = "test/bee_mocks/bold_BEEEE_2018_04_17.fa")
parser.add_argument("-L", "--libraries",		help = "paths to fastas of individual libraries/discrete samples from which zotus were found", metavar = "LIBS", default = ["test/bee_mocks/mock_"+str(i+1)+"_relabel.fa" for i in range(10)] )
parser.add_argument("-S", "--specifications",	help = "path to a text file detailing the read count binning strategy and thresholds", metavar = "SPEC", default = "specifications_test.txt")
parser.add_argument("-X", "--taxgroups",		help = "path to a two-column csv file specifying the taxon for each haplotype", metavar = "TAXA", default = "test/bee_mocks/all.dummytaxon.csv")
#parser.add_argument("-T", "--tree",			help = "path to an ultrametric tree of the zotus", metavar = "TREE", type = str, default = "test/bee_mocks/all.zotus.addglopairmaxit1k.nwk")

parser.add_argument("-o", "--output_directory",	help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-T", "--threads",			help = "number of threads to use", default = 3, metavar = "N")

	# Input file variables
parser.add_argument("-a", "--aligned",			help = "is the input fasta of unique sequences already aligned?", action = "store_true", default = True)

	# Clade finder variables
parser.add_argument("-d", "--divergence",		help = "the divergence level to use for assigning clades (default is 0.03)", default = 0.03, type = float, metavar = "N")

	# Length parameters
parser.add_argument("-n", "--minimum_length",	help = "remove sequences that are shorter than this value", type = int, default = 0, metavar = "N")
parser.add_argument("-x", "--maximum_length",	help = "remove sequences that are longer than this value", type = int, default = float('Inf'), metavar = "N")
parser.add_argument("-l", "--expected_length",	help = "the expected length of the sequences", type = int, default = 0, metavar = "N")
parser.add_argument("-p", "--percent_variation",	help = "the percentage variation from the expected length within which sequences should be kept", type = float, default = 0, metavar = "N")
parser.add_argument("-b", "--bases_variation",	help = "the number of bases of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")
parser.add_argument("-c", "--codons_variation",	help = "the number of codons of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")
parser.add_argument("--only_vary_by_codon",		help = "remove sequences that do not vary by a multiple of 3 bases from the expected length", action = "store_true")

	# Translation parameters
parser.add_argument("-s", "--table",			help = "the number referring to the translation table to use for translation filtering", metavar = "TABLE", default = 5)
parser.add_argument("-r", "--reading_frame",		help = "coding frame of sequences, if known", type = int, choices = [1,2,3], metavar = "N")
parser.add_argument("-f", "--detection_confidence",	help = "confidence level (0 < x < 1) for detection of reading frame (default 0.95, usually no need to change)", type = float, default = 0.95, metavar = "N")
parser.add_argument("-m", "--detection_minstops",	help = "minimum number of stops to encounter for detection (default 50, may need to decrease for few sequences)", type = int, default = 100, metavar = "N")

	# Reference matching parameters
parser.add_argument("-j", "--match_length",		help = "the minimum alignment length to consider a BLAST match when comparing zotus against the reference", type = int, metavar = "N", default = 350)
parser.add_argument("-i", "--match_percent",		help = "the minimum percent identity to consider a BLAST match when comparing zotus against the reference", type = int, metavar = "N", default = 97)

# Class definitions

# Function definitions

if __name__ == "__main__":
	
	#######################
	# INITIAL PREPARATION #
	#######################
	
	# Get inputs
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.zotus))[0]
	
	# Make the output directory
	
	# ?? TODO
	
	# Check options combinations
	
	#TODO: reject combination of aligned and tree - unecessary
	
	# Create dictionary of length specifications
	
	length_spec = {}
	if args.minimum_length > 0:			length_spec["min_bp"] = args.minimum_length
	if args.maximum_length < float('Inf'):	length_spec["max_bp"] = args.maximum_length
	if args.expected_length > 0:			length_spec["len_bp"] = args.expected_length
	if args.bases_variation > 0:			length_spec["var_bp"] = args.bases_variation
	if args.percent_variation > 0:			length_spec["var_pc"] = args.percent_variation
	if args.percent_variation > 0:			length_spec["var_cdn"] = args.codons_variation
	
	
	length_spec = {'len_bp' : 418,
				'var_pc' : 0}
	
	# Resolve length specifications
	
	min_max_bp = filterlength.resolve_length_spec(length_spec)
	
	###############
	# FIND CLADES #
	###############
	
	args.tree=0
	
	if(args.tree):
		tree = findclades.read_newick_string(args.tree)
		rawZotus = SeqIO.to_dict(SeqIO.parse(args.zotus, "fasta"))
		rawPath = args.zotus
	else:
		
		# Read or make the alignment
		
		if(args.aligned):
			alignedZotus = AlignIO.read(args.zotus, "fasta")
			alignedPath = args.zotus
			rawZotus = findclades.degap_alignment(alignedZotus)
			rawPath = os.path.join(args.output_directory, filename + "_unaligned.fa")
			SeqIO.write(rawZotus.values(), rawPath, "fasta")
		else:
			alignedZotus = findclades.do_alignment(args.zotus, args.threads)
			alignedPath = os.path.join(args.output_directory, filename + "_aligned.fa")
			AlignIO.write(alignedZotus, alignedPath, "fasta")
			rawZotus = findclades.degap_alignment(alignedZotus)
			rawPath = args.zotus
		
		# Compute the distance matrix and tree
		
		#dist, tree = make_dist_and_tree(align)
		tree = findclades.make_tree_R(alignedPath)
		#tree = findclades.make_tree_mafft(rawPath)
		
		del alignedZotus
		del alignedPath
	
	# Find the clades
	
	#clades = find_clades(dist, args.divergence)
	clades = findclades.get_clades_R(tree, args.divergence)
	
	# Output the tree
	
	
	if(not args.tree):
		#Phylo.write(tree, os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'newick')
		with open(os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'w') as o:
			o.write(tree)
	
	# Clean up
	
	#dist.del()
	del tree

	
	########################################
	# COMPUTE LIBRARY AND TOTAL READCOUNTS #
	########################################
	
	library_counts, total_counts = findlibraries.match_libraries_to_sequences(rawZotus, args.libraries)
	
	#############
	# READ TAXA #
	#############
	
	taxa = findtaxa.parse_taxa(args.taxgroups, rawZotus.keys())
	
	##########################
	# DESIGNATE CONTROL SETS #
	##########################
	
	# TODO: note that for now zotus refers to a dict of SeqRecords
	
	# Create length-based control list
	
	length_bad = { name for name, seq in rawZotus.items() if not filterlength.check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) }
	
	# Create translation based control list
	
		# Find reading frame
	if(args.reading_frame):
		frame = args.reading_frame
	else:
		frame = filtertranslate.detect_frame(rawZotus, args.table, args.detection_confidence, args.detection_minstops)
	
		# Filter
	trans_bad = { name for name, seq in rawZotus.items() if filtertranslate.stopcount(seq, args.table, frame) > 0 }
	
	bad = length_bad.union(trans_bad)
	
	# Create reference based control list
	
	ref_good = filterreference.blast_for_presence(rawPath, args.reference, args.match_percent, args.match_length, args.threads)
	ref_good = ref_good - bad
	
	# Finalise
	
	good = ref_good
	
	# Clean up
	
	del frame
	
	
	####################
	# CONSOLIDATE DATA #
	####################
	
	data = {"total" : total_counts,
		 "library" : library_counts,
		 "clade" : clades,
		 "taxa" : taxa}
	
	###########################################
	# READ AND PARSE FILTERING SPECIFICATIONS #
	###########################################
	
	spectext, specs, threshold_combinations = assessmentcore.parse_spec(args.specifications)
	
	###########################################
	# GENERATE COUNTS AND SCORES, THEN ASSESS #
	###########################################
	
	# Generate category counts for each specification
	
	counts = assessmentcore.counts_from_spec(specs, data)
	
	# Calculate score for threshold combination
	scores = list( assessmentcore.calc_combo_score(counts, thresh_set, good, bad, "standardised") for thresh_set in threshold_combinations )
	
	# Find threshold combination(s) matching minimum score
	
	minscore, min_thresholds, outtext = assessmentcore.get_minimum_thresholds(scores, threshold_combinations, specs)
	
	print(outtext)
	
	##################
	# OUTPUT RESULTS #
	##################
	
	# Output filtered haplotypes for threshold combination(s)
	assessmentcore.output_filtered_haplotypes(counts, min_thresholds, good, bad, rawPath, filename, args.output_directory)
	
	# Output csv of combination results for analysis
	
	
	
