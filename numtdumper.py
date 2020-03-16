#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" """

# Imports

import argparse
import os
import sys
import math
import pickle

from multiprocessing import Pool
from functools import partial

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

import categorycounting
import assessmentcore


# Global variables

parser = argparse.ArgumentParser(description = "")

	# Input / output paths
parser.add_argument("-A", "--asvs",		help = "path to a fasta of unique sequences to filter", required = True, metavar = "ASVs", type = str)
parser.add_argument("-R", "--reference",		help = "path to a fasta of known correct reference sequences", required = True, metavar = "REF", type = str)
parser.add_argument("-L", "--libraries",		help = "paths to fastas of individual libraries/discrete samples from which zotus were found", required = True , metavar = "LIB", type = str, nargs = '*')
parser.add_argument("-S", "--specification",	help = "path to a text file detailing the read count binning strategy and thresholds", required = True, metavar = "SPEC", type = str)
parser.add_argument("-X", "--taxgroups",		help = "path to a two-column csv file specifying the taxon for each haplotype", metavar = "TAXA", type = str)
parser.add_argument("-T", "--tree",		help = "path to an tree of the ASVs from a previous run", metavar = "TREE", type = str)

#Available options: g h k q u v w z

parser.add_argument("-o", "--outputdirectory",	help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-t", "--threads",		help = "number of threads to use", default = 4, metavar = "N", type = int)

	# Input file variables
parser.add_argument("-a", "--realign",		help = "force (re)alignment of the input ASVs", action = "store_true", default = False)

	# Master filtering variables
parser.add_argument("-y", "--anythreshold",	help = "reject ASVs that fail to meet any threshold (as opposed to all thresholds)", action = "store_true", default = False)
parser.add_argument("-u", "--addnull",		help = "include null thresholds to all filters (i.e. thresholds passing all reads)", action = "store_true", default = False)

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
parser.add_argument("-i", "--matchpercent",	help = "the minimum percent identity to consider a BLAST match when comparing ASVs against the reference", type = float, metavar = "N", default = 100)

# Class definitions

# Function definitions

# TODO: add more intermediate storage of results and more automated resume points

if __name__ == "__main__":
	
	#######################
	# INITIAL PREPARATION #
	#######################
	# Locate the script directory
	
	scriptdir = os.path.dirname(__file__)
	
	# Get inputs
	scriptdir = "/home/thomas/Documents/programming/bioinformatics/numtdumper/"
	os.chdir("/home/thomas/Documents/programming/bioinformatics/numtdumper_testdata/amm")
	args = parser.parse_args(['-A', '6_coleoptera_fftnsi.fasta', 
		'-R', 'dummy_reference.fasta', 
		#'-L', 'merge_fixed/10D_F_C2_.fasta', 'merge_fixed/10S_F_B5_.fasta', 'merge_fixed/11D_F_D2_.fasta', 'merge_fixed/11S_F_C5_.fasta', 'merge_fixed/12D_G_E2_.fasta', 'merge_fixed/12S_G_D5_.fasta', 'merge_fixed/13D_G_F2_.fasta', 'merge/13S_G_G6_.fasta', 'merge_fixed/14D_G_G2_.fasta', 'merge_fixed/14S_G_E5_.fasta', 'merge_fixed/15D_F_H2_.fasta', 'merge_fixed/15S_F_G5_.fasta', 'merge_fixed/16D_F_A3_.fasta', 'merge_fixed/16S_F_F5_.fasta', 'merge_fixed/17D_F_B3_.fasta', 'merge_fixed/17S_F_E6_.fasta', 'merge_fixed/18D_F_C3_.fasta', 'merge_fixed/18S_F_F6_.fasta', 'merge_fixed/19D_G_B2_.fasta', 'merge_fixed/19S_G_H5_.fasta', 'merge_fixed/1D_F_A1_.fasta', 'merge_fixed/1S_F_A4_.fasta', 'merge_fixed/20D_F_D3_.fasta', 'merge_fixed/20S_F_C6_.fasta', 'merge_fixed/21D_F_E3_.fasta', 'merge_fixed/21S_F_H6_.fasta', 'merge_fixed/22D_G_G3_.fasta', 'merge_fixed/22S_G_A6_.fasta', 'merge_fixed/23D_F_F3_.fasta', 'merge_fixed/23S_F_D6_.fasta', 'merge_fixed/24D_G_H3_.fasta', 'merge_fixed/24S_G_B6_.fasta', 'merge_fixed/2D_F_B1_.fasta', 'merge_fixed/2S_F_B4_.fasta', 'merge_fixed/3D_F_C1_.fasta', 'merge_fixed/3S_F_C4_.fasta', 'merge_fixed/4D_G_D1_.fasta', 'merge_fixed/4S_G_D4_.fasta', 'merge_fixed/5D_G_E1_.fasta', 'merge_fixed/5S_G_E4_.fasta', 'merge_fixed/6D_G_F1_.fasta', 'merge_fixed/6S_G_F4_.fasta', 'merge_fixed/7D_G_G1_.fasta', 'merge_fixed/7S_G_H4_.fasta', 'merge_fixed/8D_G_H1_.fasta', 'merge_fixed/8S_G_G4_.fasta', 'merge_fixed/9D_G_A2_.fasta', 'merge_fixed/9S_G_A5_.fasta', 'merge_fixed/N_DOM_REPS_A7_.fasta', 'merge_fixed/N_GRA_A7_.fasta', 
		'-L', 'merge/T4.fastq', 'merge/T6.fastq', 'merge/T7.fastq', 'merge/T8.fastq', 'merge/T9.fastq', 'merge/T10.fastq', 'merge/T11.fastq', 'merge/T12.fastq', 'merge/T13.fastq', 'merge/T14.fastq', 'merge/T15.fastq', 'merge/T16.fastq',
		'-S', '../../numtdumper/specifications.txt', 
		'-o', 'numtdumper/', 
		'-t', '4', 
		'-u', 
		'-l', '418', 
		'-p', '0', 
		'-s', '5', 
		'-i', '99.5'#,
	#	'-T', 'numtdumper/5_denoise_coleoptera_fftnsi_UPGMA.nwk'
		])
	
	#args = parser.parse_args()
	
	# Check inputs
	if(not all( [args.asvs, args.reference, args.libraries, args.specification] ) ):
		err = "Error, argument required to all of --asvs, --reference, --libraries and --specification"
		sys.exit(err)
	
	# Find the file name
	
	filename = os.path.splitext(os.path.basename(args.asvs))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.outputdirectory):
		os.makedirs(args.outputdirectory)
		
	
	# Check options combinations
	
	
	# Create dictionary of length specifications
	
	length_spec = {}
	if args.minimumlength > 0:			length_spec["min_bp"] = args.minimumlength
	if args.maximumlength < float('Inf'):	length_spec["max_bp"] = args.maximumlength
	if args.expectedlength > 0:			length_spec["len_bp"] = args.expectedlength
	if args.basesvariation > 0:			length_spec["var_bp"] = args.basesvariation
	if args.percentvariation > 0:		length_spec["var_pc"] = args.percentvariation
	if args.percentvariation > 0:		length_spec["var_cdn"] = args.codonsvariation
	
	# Resolve length specifications into minimum and maximum allowed
	
	min_max_bp = filterlength.resolve_length_spec(length_spec)
	
	###########################################
	# READ AND PARSE FILTERING SPECIFICATIONS #
	###########################################
	
	spectext, specs, threshold_combinations = assessmentcore.parse_spec(args.specification, args.addnull)
	
	# Check 
	
	if( not args.taxgroups and any("taxon" in e for e in spectext)):
		err = "Error, taxon specified as a binning strategy but no taxon file supplied"
		sys.exit(err)
	
	#Output specifications in table
	
	path = os.path.join(args.outputdirectory, filename + "_specifications.csv")
	with open(path, 'w') as o:
		
		# Write header
		o.write(",".join(spectext)+"\n")
		
		# Write data lines
		for thresh in threshold_combinations:
			o.write(",".join([str(t) for t in thresh])+"\n")
	
	##################
	# REPORT TO USER #
	##################
	
	sys.stdout.write("\nWelcome to NUMTdumper, let's dump those NUMTs!\n\n")
	
	sys.stdout.write("Parsed %s specifications, %s total threshold combinations\n" % (len(specs), len(threshold_combinations)) )
	
	###############
	# FIND CLADES #
	###############
	
	# Set up variables
	rawasvs = 0
	rawpath = 0
	alignedasvs = 0
	alignedpath = 0
	tree = 0
	
	# Detect alignment
	
	aligned = findclades.detect_aligned(args.asvs)
	
	# Parse in ASVs and align if necessary
	
	if(aligned and not args.realign):
		sys.stdout.write("Input ASVs detected as aligned (if this is not the case, run with the --realign option).")
		alignedasvs = AlignIO.read(args.asvs, "fasta")
		alignedpath = args.asvs
		rawasvs = findclades.degap_alignment(alignedasvs)
		rawpath = os.path.join(args.outputdirectory, filename + "_unaligned.fa")
		SeqIO.write(rawasvs.values(), rawpath, "fasta")
	else:
		sys.stdout.write("Input ASVs detected as not aligned")
		if(args.tree):
			sys.stdout.write(" but tree supplied, so not bothering to align")
			if(args.realign):
				sys.stdout.write(", --realign ignored")
			sys.stdout.write(".")
			rawasvs = SeqIO.to_dict(SeqIO.parse(args.asvs, "fasta"))
			rawpath = args.asvs
		else:
			sys.stdout.write(", running MAFFT FFT-NS-1 to align. This may take some time, skip this step by supplying an alignment to -ASVs\n")
			alignedasvs = findclades.do_alignment(args.asvs, args.threads)
			alignedpath = os.path.join(args.outputdirectory, filename + "_aligned.fa")
			AlignIO.write(alignedasvs, alignedpath, "fasta")
			rawasvs = findclades.degap_alignment(alignedpath)
			rawpath = args.asvs
	
	sys.stdout.write(" Read %s ASVs.\n" % (len(rawasvs)) )
	
	# Read in tree or build tree as required
	
	if(args.tree):
		sys.stdout.write("Reading supplied UPGMA tree from previous run.\n")
		tree = findclades.read_newick_string(args.tree)
	else:
		sys.stdout.write("Making a UPGMA tree from the alignment. This may take some time, skip this step in re-runs by supplying the tree to --tree\n")
		tree = findclades.make_tree_R(scriptdir, alignedpath, args.distancemodel)
		# Output the tree
		
		with open(os.path.join(args.outputdirectory, filename + "_UPGMA.nwk"), 'w') as o:
			o.write(tree)
	
	del alignedasvs
	del alignedpath
	
	# Find the clades
	
	sys.stdout.write("Finding clades from the tree at %s divergence.\n" % (args.divergence))
	
	clades = findclades.get_clades_R(scriptdir, tree, args.divergence)
	
	sys.stdout.write("Found %s clades from the tree at %s divergence.\n" % (len(clades), args.divergence))
	
	# Output csv of clade assignments
	
	findclades.write_clade_dict(clades, os.path.join(args.outputdirectory, filename + "_clades.csv"))
	
	# Clean up
	
	del tree
	
	#print(rawasvs)
	
	#############
	# READ TAXA #
	#############
	
	taxa = None
	
	if(args.taxgroups):
		sys.stdout.write("Reading taxa data\n")
		taxa = findtaxa.parse_taxa(args.taxgroups, rawasvs.keys())
	else:
		taxa = findtaxa.dummy_taxa(rawasvs.keys())
	
	########################################
	# COMPUTE LIBRARY AND TOTAL READCOUNTS #
	########################################
	
	sys.stdout.write("Matching library reads to ASVs to generate library ASV counts.\n")
	
	#TODO: check that all ASVs have >= 1 match (throw error if not) and warn if any libraries have no matches
	library_counts, total_counts = findlibraries.match_libraries_to_sequences(rawasvs, args.libraries)
	
	# Output csv of library counts
	categorycounting.write_count_dict(library_counts, rawasvs.keys(), os.path.join(args.outputdirectory, filename + "_ASVcounts.csv"))
	
	##########################
	# DESIGNATE CONTROL SETS #
	##########################
	
	# Create length-based control list
	
	sys.stdout.write("Identifying control non-target ASVs based on length.\n")
	
	nontarget_length = { name for name, seq in rawasvs.items() if not filterlength.check_length(seq, min_max_bp, args.expectedlength, args.onlyvarybycodon) }
	
	# Create translation based control list
	
	sys.stdout.write("Identifying control non-target ASVs based on translation.\n")
	
		# Find reading frame
	if(args.readingframe):
		frame = args.readingframe
	else:
		frame = filtertranslate.detect_frame(rawasvs, args.table, args.detectionconfidence, args.detectionminstops)
	
		# Filter
	nontarget_trans = { name for name, seq in rawasvs.items() if filtertranslate.stopcount(seq, args.table, frame) > 0 }
	
		# Finalise 
	nontarget = nontarget_length.union(nontarget_trans)
	
	if len(nontarget) > 0:
		sys.stdout.write("Found %s non-target ASVs: %s based on length variation and %s based on translation.\n" % (len(nontarget), len(nontarget_length), len(nontarget_trans)) )
	else: 
		err = "Error: no non-target ASVs could be found. Prior data filtering may have been too stringent."
		sys.exit(err)
	
	# Create reference based control list
	
	sys.stdout.write("Identifying control target ASVs based on references.\n")
		# Check if reference is alignment
	
	refaln = findclades.detect_aligned(args.reference)
	if(refaln):
		alignedref  = AlignIO.read(args.reference, "fasta")
		rawref = findclades.degap_alignment(alignedref)
		rawrefpath = os.path.join(args.outputdirectory, os.path.splitext(os.path.basename(args.reference))[0] + "_unaligned.fa")
		SeqIO.write(rawref.values(), rawrefpath, "fasta")
		args.reference = rawrefpath
	
	ref_match = filterreference.blast_for_presence(rawpath, args.outputdirectory, args.reference, args.matchpercent, args.matchlength, args.threads)
	
		# Finalise
	
	target = ref_match - nontarget
	if len(target) > 0:
		sys.stdout.write("Found %s target ASVs: %s out of all ASVs matched to reference set, %s non-target ASVs removed.\n" % (len(target), len(ref_match), len(ref_match - target)) )
	else:
		err = "Error: no target ASVs found"
		if len(ref_match) > 0:
			err = err + ", although " + len(ref_match) + " ASVs matched to reference set. Length thresholds may be too stringent and find too many non-target ASVs."
		else:
			err = err + ". Check your reference and your reference matching thresholds."
		sys.exit(err)
	
	# Clean up
	
	del frame
	
	# Output file with details of targets/non-targets
	
	
	path = os.path.join(args.outputdirectory, filename + "_control.txt")
	with open(path, 'w') as o:
		for category, asv in zip(["lengthfail"] * len(nontarget_length) + ["stopfail"] * len(nontarget_trans) + ["refpass"] * len(ref_match),
					   list(nontarget_length) + list(nontarget_trans) + list(ref_match)):
				o.write(category + "\t" + asv + "\n")
	
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
	
	sys.stdout.write("Generating binned counts\n")
	
	counts = assessmentcore.counts_from_spec(specs, data)
	counts_bak = counts
	# Calculate score for threshold combination
	
	sys.stdout.write("Assessing counts and scoring for each threshold combination.\n")
	
	stats_filename = "statsdumptemp.pydata"
	
	try:
		with open(stats_filename, "rb") as file:
			stats = pickle.load(file)
	except:
		chunksize = math.floor(len(threshold_combinations)/args.threads)
		with Pool(processes = args.threads) as pool:
			stats = pool.map(partial(assessmentcore.calc_stats, counts, set(rawasvs.keys()), target, nontarget, args.anythreshold, "standardised" ),
						 threshold_combinations, chunksize)
		with open(stats_filename, "wb") as f:
			pickle.dump(stats, f, pickle.HIGHEST_PROTOCOL)
	
	# Analyse scores to find best sets
	
	#TODO: collapse different thresholds with identical outputs.
	
	sys.stdout.write("Identifying optimal threshold sets\n")
	
	minscore, min_thresholds, outtext = assessmentcore.get_minimum_thresholds(stats, threshold_combinations, specs)
	
	#TODO: what to do if very large numbers of outputs with minimum score? Suppress outtext and output filtered ASVs for only a subset?
	
	#print(outtext)
	
	##################
	# OUTPUT RESULTS #
	##################
	
	# Output filtered haplotypes for threshold combination(s)
	
	sys.stdout.write("Writing filtered ASVs\n")
	
	#assessmentcore.output_filtered_haplotypes(counts, min_thresholds, args.anythreshold, target, nontarget, rawpath, filename, args.outputdirectory)
	
	# Output thresholds and scores
	
	sys.stdout.write("Writing all data\n")
	
	assessmentcore.write_specs_and_stats(spectext, threshold_combinations, stats, os.path.join(args.outputdirectory, filename + "_thresholds_scores_asvcounts.csv"))
	
	sys.stdout.write("\nNUMTs: dumped\n\n")
