#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by presence in reference set of sequences"""

# Imports
import os
import subprocess
import re
import sys
import pysam
import binning
from Bio import AlignIO, SeqIO


# Functions
def make_temp_bbmapwd(path, name):
    # Create temporary directory
    #TODO: create a proper temporary directory with appropriate libraries?
    outputdirectory = os.path.join(path, name)
    if not os.path.exists(outputdirectory):
        os.makedirs(outputdirectory)
    return(outputdirectory)

def get_seq_lengths(fasta_path):
    handle = open(fasta_path, 'rU')
    sequence_lengths = {}
    SeqRecords = SeqIO.parse(handle, 'fasta')
    for record in SeqRecords:   #loop through each fasta entry
        length = len(record.seq.ungap("-"))    #get sequence length
        sequence_lengths[record.id] = length
    return sequence_lengths

def refmatch_BBMap(querypath, workingdir, minlen, threads, ref_fasta, totalcounts, args,
                   fail = False):

    # lengths_query = get_seq_lengths(querypath)
    # lengths_ref = get_seq_lengths(ref_fasta)

    # Set up object to store
    BBMap_out = os.path.join(workingdir, "tempBBMap.sam")
    
    # Run BBMap
    bbmap_command = f"bbmap.sh ambig=random vslow semiperfectmode maxsites=100 ref={ref_fasta} in={querypath} out={BBMap_out} threads={threads} nodisk"
    bbmap_process = subprocess.Popen(bbmap_command, shell = True, 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.PIPE)
    bbmap_process.wait()

    # Get result
    bamfile = pysam.AlignmentFile(BBMap_out)

    # Iterate through each alignment
    dict_length_pass = {}
    for read in bamfile:
        # pick only the query sequences that were aligned
        if read.reference_name != None:
            if read.query_alignment_length >= minlen:  
                #  generate a dict to sort out query sequences that matched with the same reference sequence
                dict_length_pass[read.query_name] = read.reference_name          
    bamfile.close()

    if args.ignoreambigASVs:
        out = [v[0] for v in [[k for k in dict_length_pass if dict_length_pass[k] == v] for v in set(dict_length_pass.values())] if len(v) == 1]
    else:
        flipped = {}
        for key, value in dict_length_pass.items():
            if value not in flipped:
                flipped[value] = [key]
            else:
                flipped[value].append(key)
        
        out = []
        for key, value in flipped.items():
            if len(value) > 1:
                tmp_max = 0
                for val in value:
                    if totalcounts['total'][val] > tmp_max:
                        tmp_max = totalcounts['total'][val]
                        kept_ASV = val
                out.append(kept_ASV)
            else:
                out.append(value[0])
                 
    sys.stdout.write(f"found {len(dict_length_pass.values())} candidates\n")
    num_excluded_ASVs_nonuniq_ref = len(dict_length_pass.values()) - len(out)
    sys.stdout.write(f"Number of rejected candidates ASVs for matching the same reference sequence: {num_excluded_ASVs_nonuniq_ref}\n")
    return(set(out))
