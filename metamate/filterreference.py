#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by presence in reference set of sequences"""

# Imports
import os
import subprocess
import re
import sys
import pysam
from Bio import AlignIO, SeqIO

from . import binning

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

def get_max_mem():
    try:
        # Get total memory in bytes
        total_mem = os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES')
        # Return 85% in GB
        return int((total_mem * 0.85) / (1024**3))
    except (ValueError, AttributeError):
        return 4 # Default to 4GB if detection fails

def refmatch_BBMap(querypath, workingdir, minlen, threads, ref_fasta, totalcounts, args,
                   fail = False):

    # lengths_query = get_seq_lengths(querypath)
    # lengths_ref = get_seq_lengths(ref_fasta)

    # Set up object to store
    BBMap_out = os.path.join(workingdir, "tempBBMap.sam")
    
    # Determine mode based on expected length
    use_pacbio = False
    if args.expectedlength and args.expectedlength >= 500:
        use_pacbio = True
        
    input_fasta = querypath
    
    # Create temporary sanitized reference
    sanitized_ref = os.path.join(workingdir, "sanitized_ref.fasta")
    
    # Store mapping just in case, though we only need uniqueness
    ref_map = {} 
    
    with open(ref_fasta, "r") as handle, open(sanitized_ref, "w") as out_handle:
        for i, record in enumerate(SeqIO.parse(handle, "fasta")):
            new_id = f"Ref_{i}"
            ref_map[new_id] = record.id
            record.id = new_id
            record.description = new_id
            SeqIO.write(record, out_handle, "fasta")
            
    if use_pacbio:
        # Long Read Mode
        mem_gb = get_max_mem()
        bbmap_command = (f"mapPacBio.sh -Xmx{mem_gb}g usemodulo ambig=random semiperfectmode "
                         f"vslow maxsites=100 ref={sanitized_ref} in={input_fasta} out={BBMap_out} "
                         f"threads={threads} nodisk")
    else:
        # Standard Mode with Length Filtering
        # Filter out sequences >= 500bp
        filtered_fasta = os.path.join(workingdir, "filtered_query.fasta")
        count_kept = 0
        with open(input_fasta, "r") as handle, open(filtered_fasta, "w") as out_handle:
            for record in SeqIO.parse(handle, "fasta"):
                if len(record.seq) < 500:
                    SeqIO.write(record, out_handle, "fasta")
                    count_kept += 1
        
        if count_kept == 0:
             sys.stderr.write("Warning: No sequences < 500bp found for standard bbmap alignment.\n")
             # Create empty SAM file so that the program doesn't crash
             with open(BBMap_out, 'w') as f:
                 pass
             bbmap_command = None
        else:
            input_fasta = filtered_fasta
            bbmap_command = (f"bbmap.sh usemodulo ambig=random vslow semiperfectmode maxsites=100 "
                             f"ref={sanitized_ref} in={input_fasta} out={BBMap_out} "
                             f"threads={threads} nodisk")

    bbmap_stderr = ""
    bbmap_returncode = 0
    if bbmap_command:
        sys.stdout.write(f"Running alignment command: {bbmap_command}\n")
        bbmap_process = subprocess.Popen(bbmap_command, shell = True,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE)
        _, stderr_bytes = bbmap_process.communicate()
        bbmap_stderr = stderr_bytes.decode("utf-8", errors="replace")
        bbmap_returncode = bbmap_process.returncode

    memory_error = (
        "OutOfMemoryError" in bbmap_stderr
        or "out of memory" in bbmap_stderr.lower()
        or bbmap_returncode == 137  # SIGKILL (OOM killer)
    )

    # Get result
    dict_length_pass = {}
    # Check if SAM file is valid/non-empty before parsing
    # If we created an empty file above, or if bbmap failed, it might be empty.
    if os.path.exists(BBMap_out) and os.path.getsize(BBMap_out) > 0:
        try:
            bamfile = pysam.AlignmentFile(BBMap_out)
            # Iterate through each alignment
            for read in bamfile:
                # pick only the query sequences that were aligned
                if read.reference_name != None:
                    if read.query_alignment_length >= minlen:
                        #  generate a dict to sort out query sequences that matched with the same reference sequence
                        dict_length_pass[read.query_name] = read.reference_name
            bamfile.close()
        except ValueError:
            # Handle case where file exists but might be empty or invalid (e.g. no header)
            sys.stderr.write(f"Warning: Could not parse SAM file {BBMap_out}. Assuming no matches.\n")

    if len(dict_length_pass) == 0:
        if memory_error:
            sys.stderr.write(
                "ERROR: BBMap alignment produced no matches and a memory error was detected.\n"
                "The run was likely terminated due to insufficient memory (-Xmx limit reached).\n"
                "Try reducing the number of threads, increasing available RAM, or splitting the input.\n"
            )
        else:
            sys.stderr.write(
                "ERROR: BBMap alignment produced no matches against the reference.\n"
                "Either no query sequences matched the reference, or the alignment was terminated early due to memory limits.\n"
                "Check that the reference and query files are correct and that sufficient memory is available.\n"
            )
        sys.exit(1)

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
