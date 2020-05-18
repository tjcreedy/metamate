#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for matching library-wise sequences with master sequence list and summing reads"""

# Imports

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter
import os, sys, re

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
        for seqr in SeqIO.parse(librarypaths[0], 'fasta'):
            #seqr = next(SeqIO.parse(librarypaths[0], 'fasta'))
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
        sys.exit(f"Error: ASV(s) {', '.join(asvabsent)} were not found in the "
                  "library file or files.\n")
    
    # Check if all libraries had at least one match
    libabsent = [n for n in libnames if n not in countsbylibrary]
    if len(libabsent) > 0:
        sys.stderr.write(f"Warning: libraries {', '.join(libabsent)} had no "
                          "sequences matching any ASVs\n")
    
    # Return final dictionaries
    return countsbylibrary, countstotal
