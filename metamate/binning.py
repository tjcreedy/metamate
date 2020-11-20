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
    #TODO catch if a tsv is passed instead
    with open(taxafile, 'r') as fh:
        reader = csv.reader(fh)
        for row in reader:
            taxa[row[1]].add(row[0])
            allnames.add(row[0])
    
    # Check that all of the names are represented
    if(set(names) != allnames):
        sys.exit("Error: haplotype names in taxon file do not completely match haplotype names in zotu file")
    
    return(taxa)

def dummy_grouping(names):
    
    grouping = { 'x' : names}
    
    return(grouping)

def detect_format(path):
    with open(path) as f:
        fline = f.readline().strip()
    if len(fline)== 0:
        sys.stdout.write(f"Error: the first line of {path} is blank\n")
    
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
                libname = re.search(nameregex, seqr.id)
                if not libname:
                    sys.exit( "Error: can't detect a library name in header "
                             f"\"{seqr.id}\" on line {seqn} of "
                             f"{librarypaths[0]}\n")
                libname = libname.group(1)
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

def parse_readmap(master, mappath):
    # master, mappath = raw['asvs'], args.readmap
    
    # Generate reference set of ASV names with and without size tags
    asvnames = list(master.keys())
    asvstrip = [re.sub(';size=.*$', '', k) for k in master.keys()]
    if len(set(asvstrip)) != len(asvstrip):
        sys.stderr.write("Warning: ASV names are not unique after stripping "
                         ";size= tags, ensure your read map file uses exactly "
                         "the same ASV names as your input ASV file\n")
    
    # Set up empty dictionary for librarywise results
    countsbylibrary = dict()
    
    # Set up empty list for all zotu incidences
    asvcounts = {a: 0 for a in asvnames}
    libnames = []
    
    # Open map and determine separator and orientation
    maph = open(mappath, 'r')
    firstline = maph.readline().strip()
    sep = [s for s in [',', '\t'] if len(firstline.split(s)) > 1]
    if len(sep) == 0 or len(sep) > 1:
        sys.exit( "Error: cannot determine single unambiguous column separator"
                 f" for {mappath}, it should be either comma or tab\n")
    else:
        sep = sep[0]
    firstline = firstline.split(sep)
    if len(set(firstline)) != len(firstline):
        sys.exit(f"Error: column headings in {mappath} are not unique\n")
    asvcolumns = (any(n in firstline for n in asvnames) 
                  or any(s in firstline for s in asvstrip))
    
    # Parse if ASVs are the columns
    if asvcolumns:
        # Check that the number of ASVs matches
        if len(firstline) != len(asvnames):
            if len(firstline) == len(asvnames) + 1:
                firstline = firstline[1:]
            else:
                sys.exit(f"Error: columns of {mappath} determined to be ASVs, "
                          "but number of data columns does not mach number of "
                          "input ASVs\n")
        namessorted = []
        # Check that all ASVs are present, and if so ensure names are sorted
        if all(s in firstline for s in asvstrip):
            sindex = [asvstrip.index(f) for f in firstline]
            namessorted = [asvnames[i] for i in sindex]
        elif all(n in firstline for n in asvnames):
            namessorted = firstline
        else:
            sys.exit(f"Error: columns of {mappath} determined to be ASVs, but "
                      "cannot match all column names with all names of input "
                      "ASVs\n")
        # Parse the lines
        for n, line in enumerate(maph):
            # n, line = list(enumerate(maph))[0]
            line = line.strip().split(sep)
            if len(line) == 1:
                continue
            lib, counts = line[0], [int(c) for c in line[1:]]
            if len(counts) != len(firstline):
                sys.exit(f"Error: number of items in line {n+2} of {mappath} "
                          "does not match number of columns\n")
            if not all(type(c) is int for c in counts):
                sys.exit(f"Error: not all values of line {n+2} of {mappath} "
                          "are integers")
            countsbylibrary[lib] = dict()
            libnames.append(lib)
            for name, count in zip(namessorted, counts):
                if count > 0:
                    countsbylibrary[lib][name] = count
                    asvcounts[name] += count
    # Parse if ASVs are the rows
    else:
        libnames = firstline
        countsbylibrary = {l: dict() for l in libnames}
        readasvs = []
        for n, line in enumerate(maph):
            # n, line = list(enumerate(maph))[0]
            line = line.strip().split(sep)
            if len(line) == 1:
                continue
            name, counts = line[0], [int(c) for c in line[1:]]
            if len(counts) != len(libnames):
                if n == 0 and len(counts) + 1 == len(firstline):
                    drop, libnames = firstline[0], firstline[1:]
                    del countsbylibrary[drop]
                else:
                    sys.exit(f"Error: number of items in line {n+2} of "
                              "{mappath} does not match number of columns\n")
            if not all(type(c) is int for c in counts):
                sys.exit(f"Error: not all values of line {n+2} of {mappath} "
                          "are integers")
            if name in asvstrip:
                name = asvnames[asvstrip.index(name)]
            elif name not in asvnames:
                sys.exit(f"Error: ASV name {name} in {mappath} does not match "
                          "any input ASV names\n")
            readasvs.append(name)
            for lib, count in zip(libnames, counts):
                if count > 0:
                    countsbylibrary[lib][name] = count
                    asvcounts[name] += count
        if not all(r in asvnames for r in readasvs):
            sys.exit(f"Error: did not find count data in {mappath} for all "
                      "input ASVs")
    maph.close()
    # Check for empty libraries
    libabsent = [l for l, v in countsbylibrary.items() if len(v) == 0]
    if len(libabsent) > 0:
        sys.stderr.write(f"Warning: libraries {', '.join(libabsent)} had no "
                          "sequences matching any ASVs\n")
        countsbylibrary = {l: v for l, v in countsbylibrary.items() 
                                    if len(v) > 0}
    
    return(countsbylibrary, {'total': asvcounts})

def detect_aligned(fasta, n):
    with open(fasta) as fh:
        head, t = '', ''
        c = 0
        while c <= n:
            head += t
            t = next(fh, None)
            if t:
                if t[0] == '>': c += 1
            else:
                break
    try:
        AlignIO.read(io.StringIO(head), "fasta")
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

def make_tree_R(scriptdir, alignpath, model, threads):
    #alignpath, model, threads =  aligned['path'], args.distancemodel, args.threads
    
    # Generate script location
    scriptpath = os.path.join(scriptdir, 'maketree.R')
    
    # Run R script
    maketreecommand = subprocess.run([scriptpath, '-a', alignpath, '-m', model,
                                      '-c', str(threads)],
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
    #args, skipalign, skipmessage, outfile = [args, args.tree,", but tree supplied so need to align",                               os.path.join(args.outputdirectory, filename)]
    
    # Set up variables
    raw = dict()
    aligned = dict()
    
    # Detect alignment
    isaligned = detect_aligned(args.asvs, 1000)
    
    # Parse in ASVs and align if necessary
    if isaligned and not args.realign:
        if not skipalign:
            sys.stdout.write("Input ASVs detected as aligned (if this is not "
                             "the case, run with the --realign option). ")
        aligned['asvs'] = AlignIO.read(args.asvs, "fasta")
        aligned['path'] = args.asvs
        raw['asvs'] = degap_alignment(aligned['asvs'])
        raw['path'] = outfile + "_unaligned.fa"
        SeqIO.write(raw['asvs'].values(), raw['path'], "fasta")
    else:
        if skipalign:
            sys.stdout.write(skipmessage)
            if(args.realign):
                sys.stdout.write("Input ASVs detected as not aligned, but --"
                                 "realign ignored as alignment is unneeded. ")
            sys.stdout.flush()
            raw['asvs'] = SeqIO.to_dict(SeqIO.parse(args.asvs, "fasta"))
            raw['path'] = args.asvs
        else:
           sys.stdout.write("Input ASVs detected as not aligned, running "
                            "MAFFT FFT-NS-1 to align. This may take some "
                            "time, skip this step by supplying an alignment "
                            "to -A/--asvs\n")
           sys.stdout.flush()
           aligned['asvs'] = do_alignment(args.asvs, args.threads)
           aligned['path'] = outfile + "_aligned.fa"
           AlignIO.write(aligned['asvs'], aligned['path'], "fasta")
           raw['asvs'] = degap_alignment(aligned['asvs'])
           raw['path'] = args.asvs

    sys.stdout.write(f"Read {len(raw['asvs'])} ASVs.\n")

    return(raw, aligned)

def find_clades(args, filename):
    #filename = infilename
    
    # Locate the script directory
    scriptdir = os.path.dirname(__file__)
    
    # Get the asv dicts
    raw, aligned = parse_asvs(args, args.tree,
                              ", but tree supplied so no need to align",
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
        if len(aligned['asvs']) > 65536:
            exit("Error: the number of ASVs (over 65,536) is too high to "
                 "compute a UPGMA tree. Please run again with fewer ASVs. "
                 "We suggest removing ASVs that are highly infrequent and "
                 "fall outside of the target length or have stops in "
                 "translation. You can use the filtertranslate command to do "
                 "standalone filtering by translation\n")
        sys.stdout.flush()
        tree = make_tree_R(scriptdir, aligned['path'], args.distancemodel,
                           args.threads)
        
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