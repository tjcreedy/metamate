#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by presence in reference set of sequences"""

# Imports
import os
import subprocess

from Bio import AlignIO, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from numtdumper import binning

# Global variables


# Function definitons


def make_temp_blastwd(path, name):
    # Create temporary directory
    #TODO: create a proper temporary directory with appropriate libraries?
    outputdirectory = os.path.join(path, name)
    if not os.path.exists(outputdirectory):
        os.makedirs(outputdirectory)
    return(outputdirectory)


def make_blastdb(subjectfile, workingdir):
    
    # Check if subjectfile is aligned
    subaln = binning.detect_aligned(subjectfile)
    if(subaln):
        alnsub  = AlignIO.read(subjectfile, "fasta")
        rawsub = binning.degap_alignment(alnsub)
        name =  os.path.splitext(os.path.basename(subjectfile))[0]
        rawsubpath = os.path.join(workingdir, f"{name}_unaligned.fa")
        SeqIO.write(rawsub.values(), rawsubpath, "fasta")
        subjectfile = rawsubpath
    
    # Create blast database
    blastdb = os.path.join(workingdir, "blastdb")
    blastdbcline = f"makeblastdb -in {subjectfile} -dbtype nucl -out {blastdb}"
    blastdbprocess = subprocess.Popen(blastdbcline, shell = True, 
                                      stdout = subprocess.PIPE, 
                                      stderr = subprocess.PIPE)
    blastdbprocess.wait()
    
    return(blastdb)


def refmatch_blast(querypath, subjectdb, workingdir, percid, minlen, threads,
                   fail = False):
    #querypath, subjectdb, workingdir, percid, minlen, threads = [raw['path'], db, wd, mp, ml, args.threads]
    # Set up blast object
    blastout = os.path.join(workingdir, "tempblastout.xml")
    blastncline = NcbiblastnCommandline(query = querypath, db = subjectdb,
                                        evalue = 0.001, perc_identity = percid,
                                        num_threads = threads, outfmt = 5,
                                        out = blastout,
                                        max_target_seqs = 1000)
    #print(blastncline)
    # Run blast
    stdout, stderr = blastncline()
    
    # Get result
    blastresultfh = open(blastout)
    blastrecords = NCBIXML.parse(blastresultfh)
    
    # Work through results finding reference hits that are suitable
    
    out = []
    for blastrecord in blastrecords:
        lengthpass = []
        for alignment in blastrecord.alignments:
            for hsp in alignment.hsps:
                lengthpass.append(hsp.query_end - hsp.query_start + 1 > minlen)
        qpass = len(lengthpass) > 0 and any(lengthpass)
        if (qpass and not fail) or (not qpass and fail):
            out.append(blastrecord.query)
    
    return(set(out))
