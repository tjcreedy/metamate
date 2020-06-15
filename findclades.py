#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for assigning clades to a set of sequences"""

# Imports

from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import io
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
parser.add_argument("-T", "--threads",            help = "number of threads to use", default = 3, metavar = "N")

# Class definitions

# Function definitions

def detect_aligned(fasta):
    try:
        AlignIO.read(fasta, "fasta")
        return(True)
    except ValueError:
        return(False)

def do_alignment(fasta, threads):
    # Run MAFFT alignment
    align_cmd = MafftCommandline(input = fasta, retree = 1,
                                 maxiterate = 0, thread = int(threads))
    align_so, align_se = align_cmd()
    align = AlignIO.read(io.StringIO(align_so), "fasta")

    return (align)

#def make_dist_and_tree_py(align):
#    # Build tree
#    distance_calc = DistanceCalculator('identity')                    # Construct distance calculator object
#    dist = distance_calc.get_distance(align)                        # Build distance matrix
#    tree_builder = DistanceTreeConstructor(distance_calc, 'upgma')    # Construct tree builder object
#    tree = tree_builder.build_tree(align)                            # Build the tree
#
#    return (dist, tree)

#def find_clades_py(dist, height):
#    # Convert distance matrix to flat format
#    dist_flat = [d for row in dist.matrix[1:] for d in row[:-1]]
#
#    # Construct linkage tree
#    linkage = hierarchy.linkage(dist_flat, method = "average")
#
#    # Cut tree
#    clades = hierarchy.cut_tree(linkage, height = height)
#    clades_flat = [c for row in clades for c in row]
#
#    # Create dictionary
#    clades = defaultdict(set)
#    for clade, name in zip(clades_flat, dist.names):
#        clades[clade].add(name)
#
#
#    return(clades)

def read_newick_string(path):
    with open(path) as t:
        tree = t.read()
    return(tree)


# TODO: deal with properly addressing to the R scripts
#
#def make_tree_mafft(path):
#    # Run tree
#    subprocess.run(['mafft', '--retree', '0', '--treeout', '--averagelinkage', path])
#
#    # Get tree
#    with open(path + ".tree") as t:
#        tree = "".join([re.sub(r'^\d+_', '', l.strip()) for l in t.readlines()])
#
#    return(tree)
#
##def phylo_to_nwkstring(tree):
#

def make_tree_R(scriptdir, alignpath, model):
    # Generate script location
    scriptpath = os.path.join(scriptdir, 'maketree.R')

    # Run R script
    maketreecommand = subprocess.run([scriptpath, '-a', alignpath, '-m', model],
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE)
    # Check for errors
    mtcstderr = maketreecommand.stderr.decode("utf-8")

    if "Error" in mtcstderr:
        if('NA/NaN/Inf in foreign function call (arg 11)' in mtcstderr):
            sys.exit("Error, UPGMA tree building failed because the aligned "
                     "sequences are too different and the distance matrix "
                     "computation failed. The alignment is probably too "
                     "gappy, try making your own alignment and re-running.")
        else:
            sys.exit("Error in R UPGMA tree construction: \n" + mtcstderr)
    elif("Warning" in mtcstderr):
        sys.stderr.write(re.sub('Warning message:\n', '', mtcstderr.strip()))

    tree = maketreecommand.stdout.decode("utf-8")
    return(tree)

def get_clades_R(scriptdir, tree, height):
    # Generate script location
    getcladespath = os.path.join(scriptdir, 'getclades.R')


    if(type(tree) == str):
        tree = tree.encode()

    # Run R script
    getcladecommand = subprocess.run([getcladespath, '-i', str(height)],
                                      stdout = subprocess.PIPE,
                                      stderr = subprocess.PIPE, input = tree)
    cladestring = getcladecommand.stdout.decode("utf-8")

    # Check for errors
    gccstderr = getcladecommand.stderr.decode("utf-8")
    if "Error" in gccstderr:
        if 'the tree is not ultrametric' in gccstderr:
            sys.exit("Error, the supplied tree is not ultrametric so clades "
                     "cannot be found")
        else:
            sys.exit(f"Error in R UPGMA tree construction: {gccstderr}\n")

    # Unpack clades data
    clades = defaultdict(set)
    for line in cladestring.splitlines():
        height, name, clade = re.split('\t', line)
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
        for clade, asvs in clade_dict.items():
            for asv in asvs:
                o.write(asv + "," + clade + '\n')

def parse_asvs(args, skipalign, skipmessage, outfile):
    #args, skipalign, skipmessage, outfile = [args, args.tree,                              ", but tree supplied so need to align",                               os.path.join(args.outputdirectory, filename)]

    # Set up variables
    raw = dict()
    aligned = dict()

    # Detect alignment
    isaligned = detect_aligned(args.asvs)

    # Parse in ASVs and align if necessary
    if isaligned and not args.realign:
        sys.stdout.write("Input ASVs detected as aligned (if this is not the "
                         "case, run with the --realign option).")
        aligned['asvs'] = AlignIO.read(args.asvs, "fasta")
        aligned['path'] = args.asvs
        raw['asvs'] = degap_alignment(aligned['asvs'])
        raw['path'] = outfile + "_unaligned.fa"
        SeqIO.write(raw['asvs'].values(), raw['path'], "fasta")
    else:
        sys.stdout.write("Input ASVs detected as not aligned")
        if skipalign :
            sys.stdout.write(skipmessage)
            if(args.realign):
                sys.stdout.write(", --realign ignored")
            sys.stdout.write(".")
            raw['asvs'] = SeqIO.to_dict(SeqIO.parse(args.asvs, "fasta"))
            raw['path'] = args.asvs
        else:
           sys.stdout.write(", running MAFFT FFT-NS-1 to align. This may "
                             "take some time, skip this step by supplying an "
                             "alignment to -ASVs\n")
           aligned['asvs'] = do_alignment(args.asvs, args.threads)
           aligned['path'] = outfile + "_aligned.fa"
           AlignIO.write(aligned['asvs'], aligned['path'], "fasta")
           raw['asvs'] = degap_alignment(aligned['asvs'])
           raw['path'] = args.asvs

    sys.stdout.write(f" Read {len(raw['asvs'])} ASVs.\n")

    return(raw, aligned)

def find_clades(args, filename):

    # Locate the script directory
    scriptdir = os.path.dirname(__file__)

    # Get the asv dicts
    raw, aligned = parse_asvs(args, args.tree,
                              ", but tree supplied so need to align",
                              os.path.join(args.outputdirectory, filename))

    tree = 0

    # Read in tree or build tree as required
    if(args.tree):
        sys.stdout.write("Reading supplied UPGMA tree from previous run.\n")
        tree = read_newick_string(args.tree)
    else:
        sys.stdout.write("Making a UPGMA tree from the alignment. This may "
                         "take some time, skip this step in re-runs by "
                         "supplying the tree to --tree\n")
        tree = make_tree_R(scriptdir, aligned['path'], args.distancemodel)

        # Output the tree
        with open(os.path.join(args.outputdirectory,
                               f"{filename}_UPGMA.nwk"), 'w') as o:
            o.write(tree)

    # Find the clades

    sys.stdout.write(f"Finding clades from the tree at {args.divergence} "
                     "divergence.\n")

    clades =  get_clades_R(scriptdir, tree, args.divergence)

    sys.stdout.write(f"Found {len(clades)} clades from the tree at "
                     f"{args.divergence} divergence.\n")

    # Output csv of clade assignments

    write_clade_dict(clades, os.path.join(args.outputdirectory,
                                          f"{filename}_clades.csv"))

    return(clades, raw)




# if __name__ == "__main__":

#     # Get options

#     args = parser.parse_args()

#     # Find the file names

#     filename = os.path.splitext(os.path.basename(args.input))[0]

#     # Make the output directory

#     if not os.path.exists(args.output_directory):
#         os.makedirs(args.output_directory)

#     # Check for bad options

#     if(args.divergence > 1 or args.divergence < 0):
#         sys.ext("Error, divergence parameter should be between 0 and 1")

#     # Make the alignment

#     if(args.aligned):
#         align = AlignIO.read(args.input, "fasta")
#     else:
#         print("Aligning sequences")
#         align = do_alignment(args.input, args.threads)
#         AlignIO.write(align, os.path.join(args.output_directory, filename + "_aligned.fa"), "fasta")

#     # Compute the ditsance matrix and tree

#     print("Building tree")
#     dist, tree = make_dist_and_tree(align)

#     # Find the clades

#     print("Finding clades")
#     clades = find_clades(dist, args.divergence)

#     # Output the tree

#     Phylo.write(tree, os.path.join(args.output_directory, filename + "_UPGMA.nwk"), 'newick')

#     # Output the clades

#     csv_write = csv.writer(open(os.path.join(args.output_directory, filename + "_clades.csv"), "w"))
#     for key, val in clades.items():
#         csv_write.writerow([key, val])
#     print("Done")
