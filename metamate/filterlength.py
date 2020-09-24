#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for filtering sequences by length variants"""

# Imports

from Bio.SeqIO.FastaIO import SimpleFastaParser
import sys
import argparse
import os

# Global variables

parser = argparse.ArgumentParser(description = "Standalone tool for filtering the sequences in a multifasta according to their length. Length can be specified by an appropriate combination of the optional arguments: either 1) the expected length, 2) one or both of minimum or maximum length, or 3) a combination of the expected length and some type of variation. Filtering can be strictly limited to sequences varying by a whole number of codons from the expected length")

parser.add_argument("input", help = "input file path", metavar = "FASTA")

parser.add_argument("-o","--output_directory", help = "output directory (default is current directory)", default = "./", metavar = "OUTDIR")

parser.add_argument("-n", "--minimum_length",        help = "remove sequences that are shorter than this value", type = int, default = 0, metavar = "N")
parser.add_argument("-x", "--maximum_length",        help = "remove sequences that are longer than this value", type = int, default = float('Inf'), metavar = "N")
parser.add_argument("-l", "--expected_length",        help = "the expected length of the sequences", type = int, default = 0, metavar = "N")
parser.add_argument("-p", "--percent_variation",    help = "the percentage variation from the expected length within which sequences should be kept", type = float, default = 0, metavar = "N")
parser.add_argument("-b", "--bases_variation",        help = "the number of bases of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")
parser.add_argument("-c", "--codons_variation",        help = "the number of codons of variation from the expected length within which sequences should be kept", type = int, default = 0, metavar = "N")

parser.add_argument("-f","--onefile",            help = "rather than outputting two separate files for passing and failing sequences (the default), output sequences in one file, with the specified suffix appended to the header of failing sequences", default = False, metavar = "SUFF")
parser.add_argument("--only_vary_by_codon",        help = "remove sequences that do not vary by a multiple of 3 bases from the expected length", action = "store_true")


# Class definitons

# Function definitions

def resolve_length_spec(args, parser):
    "Resolves specifications for length filtering to just minimum and maximum"
     # Check and calculate the length specifications
    if args.expectedlength is not None:
        # If a centre given, minimum and maximum should not be given
        if args.minimumlength > 0 or args.maximumlength < float('Inf'):
            parser.error('specify either --l/--expected length, or one/both '
                         'of -n/--minimumlength, -x/--maximumlength')
        # If a centre given, find the variation specified
        nargs = 0
        bpvar = 0
        if args.basesvariation is not None:
            bpvar = args.basesvariation
            nargs += 1
        if args.percentvariation is not None:
            bpvar = round(args.percentvariation/100 * args.expectedlength)
            nargs += 1
        if args.codonsvariation is not None:
            bpvar = 3 * args.basesvariation
            nargs += 1
        nargerr = ('specify one of -b/--basesvariation, '
                   '-p/--percentvariation, -c/--codonsvariation')
        # Throw error if more than one variation is specified
        if nargs > 1:
            parser.error(f"only {nargerr}")
        # Throw error if no variation is specified but vary by codon is set
        elif nargs == 0 and args.onlyvarybycodon:
            parser.error(f"must {nargerr} if setting --onlyvarybycodon")
        # Compute the corrected minimum and maximum
        args.minimumlength = args.expectedlength - bpvar
        args.maximumlength = args.expectedlength + bpvar
    
    lset = args.minimumlength > 0 or args.maximumlength < float('Inf')
    return(args, lset)

def check_length(seqrecord, minmaxbp, lenbp, varbycodon = False):
    "Checks if a sequence record (or sequence string) passes or fails the length specifications"
    
    # Check length is ok
    if lenbp == 0 and varbycodon:
        sys.exit("Error: check_length requires lenbp > 0 if varbycodon True")
    
    # Find length
    currlen = len(seqrecord)
    
    # Determine if length falls outside specification
    out = (( minmaxbp[0] <= currlen <= minmaxbp[1] )
           and (not varbycodon or abs(lenbp - currlen) % 3 == 0))
    
    return out

def check_length_multi(seqdict, minmaxbp, args, 
                       fail = False):
    
    out = []
    for name, seqr in seqdict.items():
        seqpass = check_length(seqr, minmaxbp, args.expectedlength, 
                               args.onlyvarybycodon)
        if (seqpass and not fail) or (not seqpass and fail):
            out.append(name)
    return(out)

if __name__ == "__main__":
    
    # Get options
    
    args = parser.parse_args()
    
    # Find the file name
    
    filename = os.path.splitext(os.path.basename(args.input))[0]
    
    # Make the output directory
    
    if not os.path.exists(args.output_directory):
        os.makedirs(args.output_directory)
    
    # Check for bad options
    
    if args.only_vary_by_codon and not args.expected_length > 0:
        sys.exit("Expected length must be specified if strict codon variation is desired")
    
    min_max_bp = resolve_length_spec(args)
    
    # Output depending on options
    passcount = 0
    failcount = 0
    with open(args.input) as infasta:
        
        if not args.onefile :        # The default, two files
            
            with open(os.path.join(args.output_directory, filename + "_lengthpass.fa"), "w") as passout:
                with open(os.path.join(args.output_directory, filename + "_lengthfail.fa"), "w") as failout:
                    
                    for head, seq in SimpleFastaParser(infasta):
                        
                        if check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) :
                            passout.write(">%s\n%s\n" % (head, seq))
                            passcount += 1
                        else:
                            failout.write(">%s\n%s\n" % (head, seq))
                            failcount += 1
            
        else:
            
            with open(os.path.join(args.output_directory, filename + "_lengthfiltered.fa"), "w") as outfasta:
                
                for head, seq in SimpleFastaParser(infasta):
                    
                    if check_length(seq, min_max_bp, args.expected_length, args.only_vary_by_codon) :
                        outfasta.write(">%s\n%s\n" % (head, seq))
                        passcount += 1
                    else:
                        outfasta.write(">%s\n%s\n" % (head + args.onefile, seq))
                        failcount += 1
    
    # Print report
    print("%i of %i sequences passed length filtering" % (passcount, passcount + failcount))


