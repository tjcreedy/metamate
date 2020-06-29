#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by the presence of stop codons in translation"""

# Imports


import os
import sys
import warnings
import argparse
import scipy.stats
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

# Global variables
def getcliargs():
    parser = argparse.ArgumentParser(description = "Standalone tool for filtering the sequences in a multifasta according to whether their translation contains stop codons. All sequences must have the same reading frame. The reading frame can be supplied if known or is automatically determined. All sequences must use the same translation table, which follows the NCBI numbering convention (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)")
    
    parser.add_argument("input", metavar = "path",
                        help = "input file path")
    parser.add_argument("table", metavar = "n",
                        help = "the number referring to the translation "
                               "table to use",
                        type = int)
    parser.add_argument("-r","--reading_frame", metavar = "n",
                        help = "coding frame of sequences, if known", 
                        type = int, choices = [1,2,3])
    parser.add_argument("-o","--output_directory", metavar = "path",
                        help = "output directory (default is current "
                               " directory)", 
                        type = str, default = "./")
    parser.add_argument("-f","--onefile", metavar = "_suffix",
                        help = "rather than outputting two separate files "
                        "for passing and failing sequences (the default), "
                        "output sequences in one file, with the specified "
                        "suffix appended to the header of failing sequences", 
                        type = str, default = False)
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
    
    return(parser.parse_args())

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

def main():
    # Get options
    args = getcliargs()
    
    
    # Find the file name
    
    filename = os.path.splitext(os.path.basename(args.input))[0]
    
    # Make the output directory
    
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    
    # Check for bad options
    
    if((args.detection_confidence != 0.95 or args.detection_minstops != 100) and args.reading_frame):
        print("Warning: specifying a detection confidence or detection minstops is useless if the reading frame is known, this will be ignored")
    
    # Find frame
    
    frame = ""
    pvalue = ""
    if(args.reading_frame):
        frame = args.reading_frame
    else:
        frame = detect_frame(args.input, args.table, args.detection_confidence, args.detection_minstops)
        print("Reading frame %i detected " % (frame))
    
    
    # Output depending on options
    passcount = 0
    failcount = 0
    with open(args.input) as infasta:
        
        if not args.onefile :        # The default, two files
            
            with open(os.path.join(args.output_directory, filename + "_transpass.fa"), "w") as passout:
                with open(os.path.join(args.output_directory, filename + "_transfail.fa"), "w") as failout:
                    
                    for head, seq in SimpleFastaParser(infasta):
                        
                        if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
                            passout.write(">%s\n%s\n" % (head, seq))
                            passcount += 1
                        else:
                            failout.write(">%s\n%s\n" % (head, seq))
                            failcount += 1
            
        else:
            
            with open(os.path.join(args.output_directory, filename + "_transfiltered.fa"), "w") as outfasta:
                
                for head, seq in SimpleFastaParser(infasta):
                    
                    if stopcount(SeqRecord(Seq(seq)), args.table, frame) == 0 :
                        outfasta.write(">%s\n%s\n" % (head, seq))
                        passcount += 1
                    else:
                        outfasta.write(">%s\n%s\n" % (head + args.onefile, seq))
                        failcount += 1
    
    # Print report
    print("%i of %i sequences passed translation filtering" % (passcount, passcount + failcount))

if __name__ == "__main__":
    main()