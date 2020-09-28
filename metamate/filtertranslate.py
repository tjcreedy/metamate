#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by the presence of stop codons in translation"""

# Imports


import os
import sys
import warnings
import argparse
import scipy.stats
import copy
import textwrap as _textwrap

from Bio import SeqIO
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width,
                                                 initial_indent=indent,
                                                 subsequent_indent=indent
                                                 ) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text

class Range(argparse.Action):
    def __init__(self, minimum=None, maximum=None, *args, **kwargs):
        self.min = minimum
        self.max = maximum
        kwargs["metavar"] = "[%d-%d]" % (self.min, self.max)
        super(Range, self).__init__(*args, **kwargs)

    def __call__(self, parser, namespace, value, option_string=None):
        if not (self.min <= value <= self.max):
            msg = 'invalid choice: %r (choose from [%d-%d])' % \
                (value, self.min, self.max)
            raise argparse.ArgumentError(self, msg)
        setattr(namespace, self.dest, value)


# Function definitions

def detect_frame(seqrecords, table, pthresh = 0.95, minstops = 100):
    # Check pthresh
    if not 0 < pthresh < 1:
        sys.exit("Error in detect_frame: pthresh must be greater than 0 and"
                 "less than 1")
    # Set defaults
    counts = [0, 0, 0]
    seqrecordsrun = []
    # Check input type and assign
    if type(seqrecords) is dict:
        seqrecordsrun = iter(list(seqrecords.values()))
    elif type(seqrecords) is list:
        seqrecordsrun = iter(seqrecords)
    elif type(seqrecords) is str:
        seqrecordsrun = SeqIO.parse(seqrecords, "fasta")
    else:
        sys.exit("Error in detect_frame: type of seq_records is not dict, list"
                 "or str")
    # Do run to count stops
    p = 0
    # Check that there are at least minstops counts and the counts are
    # significantly different - if so stop, else continue
    while((sum(counts) < minstops or p < pthresh)):
        # Extract sequence object depending on input type
        try:
            seqrecord = next(seqrecordsrun)
        except:
            sys.exit("Error in detect_frame: sequences exhausted without "
                     "reaching minimum counts or significance. Try reducing "
                     "--detectionminstops or --detectionconfidence")
        # Compute counts and add to current counts itemwise
        # Count stops in all reading frames and add to current counts
        counts = [ x + y for x, y in zip(counts, stopcount(seqrecord, table))]
        # Compute p-value if enough total counts
        p = (1 - scipy.stats.chisquare(counts)[1])
    # Returns frame with minimum counts and p-value
    return counts.index(min(counts))+1

def stopcount(seqrecord, table, frame = (1 ,2, 3)):
    # Check input types
    run_frame = (frame,) if not isinstance(frame, (tuple, list)) else frame
    # Run counting
    counts = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        for i in run_frame:
            seq = seqrecord.seq[(i-1):]
            counts.append(seq.translate(table = table).count("*"))
    # Return string or list depending on length
    if(len(counts) > 1):
        return counts
    else:
        return counts[0]

def check_stops_multi(seqdict, args, fail = False):
    # Find reading frame
    if(args.readingframe):
        frame = args.readingframe
    else:
        frame = detect_frame(seqdict, args.table, args.detectionconfidence,
                             args.detectionminstops)
    
    # Filter
    out = []
    for name, seq in seqdict.items():
        stops = stopcount(seq, args.table, frame) > 0
        if (not stops and not fail) or (stops and fail):
            out.append(name)
    
    return(out)

def getcliargs(arglist = None):
    parser = argparse.ArgumentParser(
description = ("Standalone tool for filtering the sequences in a multifasta " "according to whether their translation contains stop codons. All sequences "
"must have the same reading frame relative to the start of the sequence. The "
" reading frame can be supplied if known or is automatically determined. All " "sequences must use the same translation table, which follows the NCBI "
"numbering convention "
"(https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)"))
    
    parser.add_argument("-i", "--input", metavar = "path",
                        help = "input file path", required = True)
    parser.add_argument("-t", "---table", metavar = "n",
                        help = "the number referring to the translation "
                               "table to use",
                        type = int)
    parser.add_argument("-r","--readingframe", metavar = "n",
                        help = "coding frame of sequences, if known", 
                        type = int, choices = [1,2,3])
    parser.add_argument("-o","--output", metavar = "path",
                        help = "if --outtype is 'pass', 'fail' or 'both', "
                               "the output file name, or if --outtype is "
                               "'separate', the prefix of the output file "
                               "name",
                        type = str, required = True)
    parser.add_argument("-y", "--outtype", metavar = "type",
                        help = "one of 'pass', 'fail', 'both' or 'separate', "
                               "denoting whether to output only those that "
                               "pass filtering, only those that fail "
                               "filtering, all sequences in one file adding "
                               "a';translation=pass' or ';translation=fail' "
                               "suffix in the headers, or as two separate "
                               "files with _pass or _fail suffixes to the "
                               "output file name. Default pass",
                        choices = ['pass', 'fail', 'both', 'separate'],
                        type = str, default = 'pass')
    parser.add_argument("-c","--detectionconfidence",
                        help = "confidence level for detection "
                               "of reading frame (default 0.95, usually no "
                               "need to change)",
                        type = float, default = 0.95,
                        action = Range, minimum = 0, maximum = 1)
    parser.add_argument("-m","--detectionminstops", metavar = "n", 
                        help = "minimum number of stops to encounter for "
                               "detection (default 50, may need to decrease "
                               "for few sequences)", 
                        type = int, default = 100)
    
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    
    return(args)

def main():
    # Get options
    args = getcliargs()
    #args = getcliargs(['-i', 'tests/data/3_derep.fasta', '-t', '5', '-o', 'test', '-y', 'separate'])
    # Check for bad options
    if ((args.detectionconfidence != 0.95 or args.detectionminstops != 100) 
        and args.readingframe):
        sys.stdout.write("Warning: specifying a detection confidence or "
                         "detection minstops is useless if the reading frame "
                         "is known, this will be ignored\n")
    # Load sequences
    ntseqs = list(SeqIO.parse(args.input, "fasta"))
    
    # Find frame
    frame = ""
    if(args.readingframe):
        frame = args.readingframe
    else:
        frame = detect_frame(ntseqs, args.table, args.detectionconfidence,
                             args.detectionminstops)
        sys.stdout.write(f"Reading frame {frame} detected\n")
        sys.stdout.flush()
    # Translate sequences
    aaseqs = []
    for ntr in ntseqs:
        #seqr = ntseqs[0]
        aar = copy.deepcopy(ntr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonWarning)
            aar.seq = ntr.seq[frame-1:].translate(table = args.table)
        aaseqs.append(aar)
    
    # Output depending on options
    passcount = 0
    failcount = 0
    with open(args.input) as infasta:
        if args.outtype == 'separate':        # Two files
            with open(f"{args.output}_pass.fa", "w") as passout,  open(
                    f"{args.output}_fail.fa", "w") as failout:
                for head, seq in SimpleFastaParser(infasta):
                    if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
                        passout.write(f">{head}\n{seq}\n")
                        passcount += 1
                    else:
                        failout.write(f">{head}\n{seq}\n")
                        failcount += 1
        else:
            with open(args.output, "w") as outfasta:
                for head, seq in SimpleFastaParser(infasta):
                    if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
                        head += '_pass' if args.outtype == 'both' else ''
                        if args.outtype in ['both', 'pass']:
                            outfasta.write(f">{head}\n{seq}\n")
                        passcount += 1
                    else:
                        head += '_fail' if args.outtype == 'both' else ''
                        if args.outtype in ['both', 'fail']:
                            outfasta.write(f">{head}\n{seq}\n")
                        failcount += 1
    
    # Print report
    sys.stdout.write(f"{passcount} of {passcount + failcount} sequences "
                     "passed translation filtering\n")

if __name__ == "__main__":
    main()