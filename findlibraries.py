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
def count_asvs_in_libraries(master, librarypaths):
    # master, library_paths = [raw['asvs'], args.libraries]
    "Work through libraries counting incidences of each master sequence"
    
    # Convert master into a simple string dictionary
    asvseqs = {str(r.seq).upper() : n for n, r in master.items()}
    
    # TODO: ensure number of items in master_dict matches number of items in master
    
    # Set up empty dictionary for librarywise results
    countsbylibrary = dict()
    
    # Set up empty list for all zotu incidences
    asvcounts = []
    
    # Loop through libraries
    for path in librarypaths:
        #path = librarypaths[0]
        
        # Extract library name
        libname = os.path.splitext(os.path.basename(path))[0]
        
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
                exit()
        
        # Get list of all zotus in library
        libasvs = [asvseqs[seq] for seq in seqs if seq in asvseqs]
        
        # Add to librarywise results
        countsbylibrary[libname] = Counter(libasvs)
        
        # Add to master zotu incidences
        asvcounts.extend(libasvs)
    
    # Set up dictionary for totals results
    countstotal = {'total': Counter(asvcounts)}
    
    # Return final dictionaries
    return countsbylibrary, countstotal
