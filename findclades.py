#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for assigning clades to a set of sequences"""

# Imports

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from io import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.cluster import hierarchy
import os
import argparse
import csv
import sys
import subprocess
import re
import warnings
from collections import defaultdict

#from collections import defaultdict

# Global variables
parser = argparse.ArgumentParser(description = "Standalone tool for assigning a set of sequences to clades according to a specified divergence. Input is a set of sequences, aligned or not. If an unaligned set of sequences is supplied, sequences are aligned using MAFFT. A tree is built using UPGMA, and the tree and a text file recording clade assignments for each sequence are output to the specified directory")

parser.add_argument("input", help = "input file path", metavar = "FASTA")
parser.add_argument("-o", "--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")
parser.add_argument("-a", "--aligned", help = "is the input fasta already aligned?", action = "store_true")
parser.add_argument("-d", "--divergence", help = "the divergence level to use for assigning clades (default is 0.03)", default = 0.03, type = float, metavar = "N")
parser.add_argument("-T", "--threads",			help = "number of threads to use", default = 3, metavar = "N")

# Class definitions

# Function definitions

def do_alignment(fasta, threads):
	# Run MAFFT alignment
	align_cmd = MafftCommandline(input = fasta, retree = 2, maxiterate = 2, threads = threads)	# Construct MAFFT commandline object
	align_so, align_se = align_cmd()															# Run the commandline object
	align = AlignIO.read(StringIO(align_so), "fasta")											# Parse the standardout into alignment object
	
	return (align)

def make_dist_and_tree_py(align):
	# Build tree
	distance_calc = DistanceCalculator('identity')					# Construct distance calculator object
	dist = distance_calc.get_distance(align)						# Build distance matrix
	tree_builder = DistanceTreeConstructor(distance_calc, 'upgma')	# Construct tree builder object
	tree = tree_builder.build_tree(align)							# Build the tree
	
	return (dist, tree)

def find_clades_py(dist, height):
	# Convert distance matrix to flat format
	dist_flat = [d for row in dist.matrix[1:] for d in row[:-1]]
	
	# Construct linkage tree
	linkage = hierarchy.linkage(dist_flat, method = "average")
	
	# Cut tree
	clades = hierarchy.cut_tree(linkage, height = height)
	clades_flat = [c for row in clades for c in row]
	
	# Create dictionary
	clades = defaultdict(set)
	for clade, name in zip(clades_flat, dist.names):
		clades[clade].add(name)
	
	
	return(clades)

def read_newick_string(path):
	with open(path) as t:
		tree = t.read()
	return(tree)


# TODO: deal with properly addressing to the R scripts

def make_tree_mafft(path):
	# Run tree
	subprocess.run(['mafft', '--retree', '0', '--treeout', '--averagelinkage', path])
	
	# Get tree
	with open(path + ".tree") as t:
		tree = "".join([re.sub(r'^\d+_', '', l.strip()) for l in t.readlines()])
	
	return(tree)

#def phylo_to_nwkstring(tree):
	

def make_tree_R(alignpath, model):
	# Run R script
	maketreecommand = subprocess.run(['./maketree.R', '-a', alignpath, '-m', model], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	
	# Check for errors
	mtcstderr = maketreecommand.stderr.decode("utf-8")
	
	if("Error" in mtcstderr):
		if('NA/NaN/Inf in foreign function call (arg 11)' in mtcstderr):
			sys.exit("Error, UPGMA tree building failed because the aligned sequences are too different and the distance matrix computation failed. The alignment is probably too gappy, try making your own alignment and re-running.")
		else:
			sys.exit("Error in R UPGMA tree construction: \n" + mtcstderr)
	elif("Warning" in mtcstderr):
		warnings.warn(re.sub(r'Warning message:\n', '', mtcstderr.strip()))
	
	
	tree = maketreecommand.stdout.decode("utf-8")
	return(tree)

def get_clades_R(tree, height):
	if(type(tree) == str):
		tree = tree.encode()
	
	# Run R script
	getcladecommand = subprocess.run(['./getclades.R', '-i', str(height)], stdout = subprocess.PIPE, stderr = subprocess.PIPE, input = tree)
	cladestring = getcladecommand.stdout.decode("utf-8")
	
	if("Error" in getcladecommand.stderr.decode("utf-8")):
		if('the tree is not ultrametric' in getcladecommand.stderr.decode("utf-8")):
			sys.exit("Error, the supplied tree is not ultrametric so clades cannot be found")
		else:
			sys.exit("Error in R UPGMA tree construction: \n" + getcladecommand.stderr.decode("utf-8"))
	
	# Unpack clades data
	clades = defaultdict(set)
	for line in cladestring.splitlines():
		height, name, clade = re.split(r'\t', line)
		clades[clade].add(name)
	
	return(clades)

def degap_alignment(alignment):
	degapped = []
	for record in alignment:
		degap_record = record
		degap_record.seq = Seq(re.sub('-', '', str(degap_record.seq)).upper())
		degapped.append(degap_record)
	
	return(SeqIO.to_dict(degapped))

def write_clade_dict(clade_dict, path):
	with open(path, "w") as o:
		for clade, asv in clades.items():
			o.write([key, val,'\n'])


if __name__ == "__main__":
	
	# Get options
	
	args = parser.parse_args()
	
	# Find the file names
	
	filename = os.path.splitext(os.path.basename(args.input))[0]
	
	# Make the output directory
	
	if not os.path.exists(args.output_directory):
		os.makedirs(args.output_directory)
	
	# Check for bad options
	
	if(args.divergence > 1 or args.divergence < 0):
		sys.ext("Error, divergence parameter should be between 0 and 1")
	
	# Make the alignment
	
	if(args.aligned):
		align = AlignIO.read(args.input, "fasta")
	else:
		print("Aligning sequences")
		align = do_alignment(args.input, args.threads)
		AlignIO.write(align, os.path.join(args.output_directory, filename + "_aligned.fa"), "fasta")
	
	# Compute the ditsance matrix and tree
	
	print("Building tree")
	dist, tree = make_dist_and_tree(align)
	
	# Find the clades
	
	print("Finding clades")
	clades = find_clades(dist, args.divergence)
	
	# Output the tree
	
	Phylo.write(tree, os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'newick')
	
	# Output the clades
	
	csv_write = csv.writer(open(os.path.join(args.output_directory, filename + "_clades.csv"), "w"))
	for key, val in clades.items():
		csv_write.writerow([key, val])
	print("Done")
