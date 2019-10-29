#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" """

# Imports

#import sys
import argparse
import os

from Bio import AlignIO
#from Bio import Phylo
from Bio import SeqIO
#from Bio.SeqIO.FastaIO import SimpleFastaParser

#import re

import filterlength
import filterreference
import filtertranslate

import findlibraries
import findclades
import findtaxa

#import categorycounting
import assessmentcore


# Global variables

parser = argparse.ArgumentParser(description = "")

	# Input / output paths
parser.add_argument("-A", "--ASVs",		help = "path to a fasta of unique sequences to filter", metavar = "ASVs", type = str, default = "test/basetest/zotus.aln.fa")
parser.add_argument("-R", "--reference",		help = "path to a fasta of known correct reference sequences", metavar = "REF", type = str, default = "test/basetest/bold_BEEEE_2018_04_17.fa")
parser.add_argument("-L", "--libraries",		help = "paths to fastas of individual libraries/discrete samples from which zotus were found", metavar = "LIBS", default = ["test/basetest/sample"+str(i+1)+".fa" for i in range(10)] )
parser.add_argument("-S", "--specifications",	help = "path to a text file detailing the read count binning strategy and thresholds", metavar = "SPEC", default = "specifications_test.txt")
parser.add_argument("-X", "--taxgroups",		help = "path to a two-column csv file specifying the taxon for each haplotype", metavar = "TAXA", default = "test/basetest/zotus.dummytaxon.csv")
#parser.add_argument("-T", "--tree",			help = "path to an ultrametric tree of the zotus", metavar = "TREE", type = str, default = "test/bee_mocks/all.zotus.addglopairmaxit1k.nwk")

Available options: g h k q u v w z

parser.add_argument("-o", "--outputdirectory",	help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-t", "--threads",		help = "number of threads to use", default = 4, metavar = "N")

	# Input file variables
parser.add_argument("-a", "--aligned",		help = "is the input fasta of unique sequences already aligned?", action = "store_true", default = True)

	# Master filtering variables
parser.add_argument("-y", "--anythreshold",	help = "reject ASVs that fail to meet any threshold (as opposed to all thresholds)", action = "store_true", default = False)
parser.add_argument("-u", "--addnull",		help = "include null thresholds to all filters (i.e. thresholds passing all reads)", action = "stor_true", default = False)

	# Clade finder variables
parser.add_argument("-m", "--distancemodel",	help = "substitution model for UPGMA tree estimation (passed to R dist.dna, default F84)", default = "F84", type = str, metavar = "X")
parser.add_argument("-d", "--divergence",	help = "divergence level to use for assigning clades (default is 0.2)", default = 0.2, type = float, metavar = "N")

	# Length parameters
parser.add_argument("-n", "--minimumlength",	help = "designate ASVs that are shorter than this value as non-target", type = int, default = 0, metavar = "N")
parser.add_argument("-x", "--maximumlength",	help = "designate ASVs that are longer than this value as non-target", type = int, default = float('Inf'), metavar = "N")
parser.add_argument("-l", "--expectedlength",	help = "the expected length of the sequences", type = int, default = 0, metavar = "N")
parser.add_argument("-p", "--percentvariation",	help = "the percentage variation from the expected length within which ASVs should not be designated as non-target", type = float, default = 0, metavar = "N")
parser.add_argument("-b", "--basesvariation",	help = "the number of bases of variation from the expected length within which ASVs should not be designated as non-target", type = int, default = 0, metavar = "N")
parser.add_argument("-c", "--codonsvariation",	help = "the number of codons of variation from the expected length within which ASVs should not be designated as non-target", type = int, default = 0, metavar = "N")
parser.add_argument("--onlyvarybycodon",		help = "designate ASVs that do not vary by a multiple of 3 bases from the expected length as non-target", action = "store_true")

	# Translation parameters
parser.add_argument("-s", "--table",		help = "the number referring to the translation table to use for translation filtering", metavar = "TABLE", default = 5)
parser.add_argument("-r", "--readingframe",	help = "coding frame of sequences, if known", type = int, choices = [1,2,3], metavar = "N")
parser.add_argument("-f", "--detectionconfidence",	help = "confidence level (0 < x < 1) for detection of reading frame (default 0.95, usually no need to change)", type = float, default = 0.95, metavar = "N")
parser.add_argument("-e", "--detectionminstops",	help = "minimum number of stops to encounter for detection (default 100, may need to decrease for few sequences)", type = int, default = 100, metavar = "N")

	# Reference matching parameters
parser.add_argument("-j", "--matchlength",	help = "the minimum alignment length to consider a BLAST match when comparing ASVs against the reference", type = int, metavar = "N", default = 350)
parser.add_argument("-i", "--matchpercent",	help = "the minimum percent identity to consider a BLAST match when comparing ASVs against the reference", type = int, metavar = "N", default = 97)

# Class definitions

# Function definitions

if __name__ == "__main__":
	
	#######################
	# INITIAL PREPARATION #
	#######################
	
	# Get inputs
	
	args = parser.parse_args()
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.asvs))[0]
	
	# Make the output directory
	
	# ?? TODO
	
	# Check options combinations
	
	#TODO: reject combination of aligned and tree - unecessary
	
	# Create dictionary of length specifications
	
	length_spec = {}
	if args.minimumlength > 0:		length_spec["min_bp"] = args.minimumlength
	if args.maximumlength < float('Inf'):	length_spec["max_bp"] = args.maximumlength
	if args.expectedlength > 0:		length_spec["len_bp"] = args.expectedlength
	if args.basesvariation > 0:		length_spec["var_bp"] = args.basesvariation
	if args.percentvariation > 0:		length_spec["var_pc"] = args.percentvariation
	if args.percentvariation > 0:		length_spec["var_cdn"] = args.codonsvariation
	
	
	length_spec = {'len_bp' : 418,
			  'var_pc' : 0}
	
	# Resolve length specifications into minimum and maximum allowed
	
	min_max_bp = filterlength.resolve_length_spec(length_spec)
	
	###########################################
	# READ AND PARSE FILTERING SPECIFICATIONS #
	###########################################
	
	spectext, specs, threshold_combinations = assessmentcore.parse_spec(args.specifications, args.addnull)
	
	#############
	# READ TAXA #
	#############
	
	taxa = None
	
	if(args.taxgroups):
		taxa = findtaxa.parse_taxa(args.taxgroups, rawasvs.keys())
	else:
		if(any("taxon" in e for e in spectext)):
			err = "Error, taxon specified as a binning strategy but no taxon file supplied"
			sys.exit(err)
		else:
			taxa = findtaxa.dummy_taxa(rawasvs.keys())
	
	###############
	# FIND CLADES #
	###############
	
	args.tree=0
	
	if(args.tree):
		tree = findclades.read_newick_string(args.tree)
		rawasvs = SeqIO.to_dict(SeqIO.parse(args.asvs, "fasta"))
		rawpath = args.asvs
	else:
		
		# Read or make the alignment
		
		if(args.aligned):
			alignedasvs = AlignIO.read(args.asvs, "fasta")
			alignedpath = args.asvs
			rawasvs = findclades.degap_alignment(alignedasvs)
			rawpath = os.path.join(args.output_directory, filename + "_unaligned.fa")
			SeqIO.write(rawasvs.values(), rawpath, "fasta")
		else:
			alignedasvs = findclades.do_alignment(args.asvs, args.threads)
			alignedpath = os.path.join(args.output_directory, filename + "_aligned.fa")
			AlignIO.write(alignedasvs, alignedpath, "fasta")
			rawasvs = findclades.degap_alignment(alignedpath)
			rawasvs = args.asvs
		
		# Compute the distance matrix and tree
		
		#dist, tree = make_dist_and_tree(align)
		tree = findclades.make_tree_R(alignedpath, args.model)
		#tree = findclades.make_tree_mafft(rawPath)
		
		del alignedasvs
		del alignedpath
	
	# Find the clades
	
	#clades = find_clades(dist, args.divergence)
	clades = findclades.get_clades_R(tree, args.divergence)
	
	# Output the tree
	
	if(not args.tree):
		#Phylo.write(tree, os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'newick')
		with open(os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'w') as o:
			o.write(tree)
	
	# Output csv of clade assignments
	
	findclades.write_clade_dict(clades, os.path.join(args.output_directory, filename + "_clades.csv"))
	
	# Clean up
	
	#dist.del()
	del tree
	
	########################################
	# COMPUTE LIBRARY AND TOTAL READCOUNTS #
	########################################
	
	library_counts, total_counts = findlibraries.match_libraries_to_sequences(rawasvs, args.libraries)
	
	# Output csv of library counts
	categorycounting.write_count_dict(library_counts, rawasvs.keys(), os.path.join(args.output_directory, filename + "_ASVcounts.csv"))
	
	##########################
	# DESIGNATE CONTROL SETS #
	##########################
	
	# Create length-based control list
	
	nontarget_length = { name for name, seq in rawasvs.items() if not filterlength.check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) }
	
	# Create translation based control list
	
		# Find reading frame
	if(args.reading_frame):
		frame = args.reading_frame
	else:
		frame = filtertranslate.detect_frame(rawasvs, args.table, args.detection_confidence, args.detection_minstops)
	
		# Filter
	nontarget_trans = { name for name, seq in rawasvs.items() if filtertranslate.stopcount(seq, args.table, frame) > 0 }
	
	nontarget = nontarget_length.union(nontarget_trans)
	
	# Create reference based control list
	
	target_ref = filterreference.blast_for_presence(rawpath, args.reference, args.match_percent, args.match_length, args.threads)
	target_ref = target_ref - nontarget
	
	# Finalise
	
	target = target_ref
	
	# Clean up
	
	del frame
	
	####################
	# CONSOLIDATE DATA #
	####################
	
	data = {"total" : total_counts,
		 "library" : library_counts,
		 "clade" : clades,
		 "taxon" : taxa}
	
	###########################################
	# GENERATE COUNTS AND SCORES, THEN ASSESS #
	###########################################
	
	# Generate category counts for each specification
	
	counts = assessmentcore.counts_from_spec(specs, data)
	
	# Calculate score for threshold combination
	
	scores = list( assessmentcore.calc_stats(counts, threshold_set, rawasvs.keys(), target, nontarget, args.anythreshold, "standardised") for threshold_set in threshold_combinations )
	
	# Find threshold combination(s) matching minimum score
	
	minscore, min_thresholds, outtext = assessmentcore.get_minimum_thresholds(scores, threshold_combinations, specs)
	
	print(outtext)
	
	##################
	# OUTPUT RESULTS #
	##################
	
	# Output filtered haplotypes for threshold combination(s)
	assessmentcore.output_filtered_haplotypes(counts, min_thresholds, target, nontarget, rawpath, filename, args.output_directory)
	
	# Output thresholds and scores
	assessmentcore.write_specs_and_stats(spectext, threshold_combinations, scores, os.path.join(args.output_directory, filename + "_thresholds_scores_asvcounts.nwk"))

