#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for assigning clades to a set of sequences"""

# Imports


import os
import sys
import subprocess
import re
import io
import csv

from collections import defaultdict, Counter
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.Align.Applications import MafftCommandline

# Function definitions


def parse_taxa(taxafile, names):
    
    taxa = defaultdict(set)
    allnames = set()
    
    # Load in taxon file csv to dict
    with open(taxafile, 'r') as fh:
        reader = csv.reader(fh)
        for row in reader:
            taxa[row[1]].add(row[0])
            allnames.add(row[0])
    
    # Check that all of the names are represented
    if(set(names) != allnames):
        sys.exit("Error: haplotype names in taxon file do not completely match haplotype names in zotu file")
    
    return(taxa)

def dummy_taxa(names):
    
    taxa = { 'x' : names}
    
    return(taxa)

def detect_format(path):
    with open(path) as f:
        fline = f.readline().strip()
    
    if(fline[0] == "@"):
        return("fastq")
    elif(fline[0] == ">"):
        return("fasta")
    else:
        return("unknown")

def count_asvs_in_libraries(master, librarypaths):
    # master, librarypaths = [raw['asvs'], args.libraries]
    "Work through libraries counting incidences of each master sequence"
    
    # Convert master into a simple string dictionary
    asvseqs = {str(r.seq).upper() : n for n, r in master.items()}
    
    # TODO: ensure number of items in master_dict matches number of items in master
    
    # Set up empty dictionary for librarywise results
    countsbylibrary = dict()
    
    # Set up empty list for all zotu incidences
    asvcounts = []
    
    libnames = []
    
    if len(librarypaths) > 1:
        # Loop through libraries
        for path in librarypaths:
            #path = librarypaths[0]
            # Extract library name
            libname = os.path.splitext(os.path.basename(path))[0]
            libnames.append(libname)
            # Get sequences
            seqformat = detect_format(path)
            seqs = []
            with open(path) as fh:
                if seqformat == 'fasta':
                    seqs = [s.upper() for h, s in SimpleFastaParser(fh)]
                elif seqformat == 'fastq':
                    seqs = [s.upper() for h, s, q in FastqGeneralIterator(fh)]
                else:
                    sys.exit(f"Error: can't detect format of {path}\n")
            # Get list of all zotus in library
            libasvs = [asvseqs[seq] for seq in seqs if seq in asvseqs]
            # Add to librarywise results
            countsbylibrary[libname] = Counter(libasvs)
            # Add to master zotu incidences
            asvcounts.extend(libasvs)
    else:
        nameregex = ";(?:barcodelabel|sample)=([^;]+);"
        seqn = 0
        seqformat = detect_format(librarypaths[0])
        for seqr in SeqIO.parse(librarypaths[0], seqformat):
            #seqr = next(SeqIO.parse(librarypaths[0], seqformat))
            seqn += 1
            seq = str(seqr.seq).upper()
            if seq in asvseqs:
                libname = re.search(nameregex, seqr.id).group(1)
                if not libname:
                    sys.exit( "Error: can't detect a library name in header "
                             f"\"{seqr.id}\" on line {seqn} of "
                             f"{librarypaths[0]}\n")
                asvname = asvseqs[seq]
                if libname not in countsbylibrary:
                    countsbylibrary[libname] = {asvname: 1}
                elif asvname not in countsbylibrary[libname]:
                    countsbylibrary[libname][asvname] = 1
                else:
                    countsbylibrary[libname][asvname] += 1
                asvcounts.append(asvname)
    
    # Set up dictionary for totals results
    countstotal = {'total': Counter(asvcounts)}
    
    # Check if all ASVs were in at least one master
    asvabsent = [n for n, r in master.items() if n not in countstotal['total']]
    if len(asvabsent) > 0:
        sys.exit(f"Error: {len(asvabsent)} ASV(s) were not found in the "
                 f"library file(s):\n{', '.join(asvabsent)}\n")
    
    # Check if all libraries had at least one match
    libabsent = [n for n in libnames if n not in countsbylibrary]
    if len(libabsent) > 0:
        sys.stderr.write(f"Warning: libraries {', '.join(libabsent)} had no "
                          "sequences matching any ASVs\n")
    
    # Return final dictionaries
    return countsbylibrary, countstotal



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

def read_newick_string(path):
    with open(path) as t:
        tree = t.read()
    return(tree)


# TODO: deal with properly addressing to the R scripts

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
                             "alignment to -A/--asvs\n")
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
                         "supplying the tree to -T/--tree\n")
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