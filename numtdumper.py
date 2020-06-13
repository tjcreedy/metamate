#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" """

# Imports

import argparse
import os
import sys
import math
import textwrap as _textwrap
from multiprocessing import Pool
from functools import partial

import filterlength
import findlibraries
import findclades
import findtaxa

import categorycounting
import assessmentcore

# Class definitions

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

# Function definitions


# Global variables


# B D E F H I J K M N O P Q U V W Y Z
# e f g h i j k m q u v w z

parser = argparse.ArgumentParser(description = "")
parser._positionals.title = "action to perform"
parser._optionals.title = "optional arguments"


parser.add_argument("-t", "--threads",
                        help = "number of threads to use (default 1)", 
                        default = 1, metavar = "N", type = int)

subparsers = parser.add_subparsers(help = '', dest = 'action')

    # Core inputs on a separate parser which acts as parent to subs

coreparser = argparse.ArgumentParser(add_help = False)
core = coreparser.add_argument_group('core inputs and options')
core.add_argument("-A", "--asvs",
                  help = "path to a fasta of unique sequences to filter",
                  required = True, metavar = "ASVs", type = str)
core.add_argument("-L", "--libraries",
                  help = "path to fastx file(s) of individual libraries"
                         "/discrete samples from which ASVs were drawn, or "
                         "a single fastx with ;samplename=.*; or "
                         ";barcodelabel=.*; annotations in headers.",
                  metavar = "LIB", type = str, nargs = '*')
core.add_argument("-G", "--taxgroups",
                  help = "path to a two-column csv file specifying the "
                         "taxon for each ASV", 
                  metavar = "TAXA", type = str)
core.add_argument("-T", "--tree",
                  help = "path to an tree of the ASVs from a previous run",
                  metavar = "TREE", type = str)
core.add_argument("-o", "--outputdirectory",
                  help = "output directory (default is current directory)",
                  default = "./", metavar = "OUTDIR")
core.add_argument("-a", "--realign",
                  help = "force (re)alignment of the input ASVs", 
                  action = "store_true", default = False)
core.add_argument("-y", "--anyfail",
                  help = "reject ASVs when any incidences fail to meet "
                         "a threshold (default is all incidences)", 
                  action = "store_true", default = False)
#TODO: throw error if not A + C or A + R + L + S 

findparser = subparsers.add_parser("find", parents = [coreparser],
                                   help = "find NUMTs by supplying threshold "
                                          "ranges and control specifications")
findparser._optionals.title = "arguments"

findparser.add_argument("-S", "--specification",
                        help = "path to a text file detailing the read count "
                               "binning strategy and thresholds",
                        required = True, metavar = "SPEC", type = str)

#parser.add_argument("-u", "--addnull",
#                    help = "include null thresholds to all filters (i.e. "
#                    "thresholds passing all reads)", 
#                    action = "store_true", default = False)

    # Clade finder variables
cladfind = findparser.add_argument_group('tlade finding')
cladfind.add_argument("--distancemodel",
                      help = "substitution model for UPGMA tree estimation "
                             "(passed to R dist.dna, default F84)", 
                      default = "F84", type = str, metavar = "X")
cladfind.add_argument("-d", "--divergence",
                      help = "divergence level to use for assigning clades "
                             "(default is 0.2)",
                      default = 0.2, type = float, metavar = "N")

    # Reference matching parameters
refmatch = findparser.add_argument_group("teference-matching-based target "
                                         "identification")
refmatch.add_argument("-R", "--references",
                      help = "path to a fasta of known correct reference "
                             "sequences", 
                      metavar = "REF", type = str)
refmatch.add_argument("--matchlength",
                      help = "the minimum alignment length to consider a "
                             "BLAST match when comparing ASVs against "
                             "reference sequences",
                      type = int, metavar = "N", default = 350)
refmatch.add_argument("--matchpercent",
                      help = "the minimum percent identity to consider a "
                             "BLAST match when comparing ASVs against "
                             "reference sequences", 
                      type = float, metavar = "N", default = 100)

    # Length parameters
lengths = findparser.add_argument_group("tength-based non-target "
                                        "identification")
lengths.add_argument("-n", "--minimumlength",
                     help = "designate ASVs that are shorter than this value "
                            "as non-target", 
                     type = int, default = 0, metavar = "N")
lengths.add_argument("-x", "--maximumlength",
                     help = "designate ASVs that are longer than this value "
                            "as non-target",
                     type = int, default = float('Inf'), metavar = "N")
lengths.add_argument("-l", "--expectedlength",
                     help = "the expected length of the sequences", 
                     type = int, default = 0, metavar = "N")
lengths.add_argument("-p", "--percentvariation",
                     help = "the percentage variation from the expected "
                            "length within which ASVs should not be "
                            "designated as non-target", 
                     type = float, default = 0, metavar = "N")
lengths.add_argument("-b", "--basesvariation",
                     help = "the number of bases of variation from the "
                            "expected length within which ASVs should not be "
                            "designated as non-target", 
                    type = int, default = 0, metavar = "N")
lengths.add_argument("-c", "--codonsvariation",
                     help = "the number of codons of variation from the "
                            "expected length within which ASVs should not be "
                            "designated as non-target", 
                     type = int, default = 0, metavar = "N")
lengths.add_argument("--onlyvarybycodon",
                     help = "designate ASVs that do not vary by a multiple of "
                            "3 bases from the expected length as non-target",
                     action = "store_true")

    # Translation parameters
transl = findparser.add_argument_group("translation-based non-target "
                                       "identification")
transl.add_argument("-s", "--table",
                    help = "the number referring to the translation table "
                           "to use for translation filtering", 
                    metavar = "TABLE", default = 5)
transl.add_argument("-r", "--readingframe",
                    help = "coding frame of sequences, if known", 
                    type = int, choices = [1,2,3], metavar = "N")
transl.add_argument("--detectionconfidence",
                    help = "confidence level (0 < x < 1) for detection of "
                           "reading frame (default 0.95, usually no need to "
                           "change)",
                    type = float, default = 0.95, metavar = "N")
transl.add_argument("--detectionminstops",
                    help = "minimum number of stops to encounter for "
                           "detection (default 100, may need to decrease for "
                           "few input ASVs)",
                    type = int, default = 100, metavar = "N")

dumpparser = subparsers.add_parser("dump", parents = [coreparser],
                                   help = "dump NUMTs found in a previous run "
                                          "or with fixed thresholds")
dumpparser._optionals.title = "arguments"
#TODO: make sure resultcache.json is correct name
dumpparser.add_argument("-C", "--resultcache",
                        help = "path to the resultcache.json file from a "
                               "previous run ",
                        metavar = 'CACHE', type = str)
dumpparser.add_argument("-S", "--specification",
                        help = "one or more [category(/ies); metric; "
                               "threshold] strings denoting the specification "
                               "for dumping NUMTs",
                        metavar = "[C(s); M; T]", type = str, nargs = '*')





# Class definitions

# Function definitions

# TODO: add more intermediate storage of results and more automated resume points

if __name__ == "__main__":
    
    #######################
    # INITIAL PREPARATION #
    #######################
    
    
    # Get inputs
#    scriptdir = "/home/thomas/Documents/programming/bioinformatics/numtdumper/"
#    os.chdir("/home/thomas/seqtesting/NUMTdumper/amm")
#    args = parser.parse_args(['-A', '6_coleoptera_fftnsi.fasta', 
#        '-R', 'dummy_reference.fasta', 
#        #'-L', 'merge_fixed/10D_F_C2_.fasta', 'merge_fixed/10S_F_B5_.fasta', 'merge_fixed/11D_F_D2_.fasta', 'merge_fixed/11S_F_C5_.fasta', 'merge_fixed/12D_G_E2_.fasta', 'merge_fixed/12S_G_D5_.fasta', 'merge_fixed/13D_G_F2_.fasta', 'merge/13S_G_G6_.fasta', 'merge_fixed/14D_G_G2_.fasta', 'merge_fixed/14S_G_E5_.fasta', 'merge_fixed/15D_F_H2_.fasta', 'merge_fixed/15S_F_G5_.fasta', 'merge_fixed/16D_F_A3_.fasta', 'merge_fixed/16S_F_F5_.fasta', 'merge_fixed/17D_F_B3_.fasta', 'merge_fixed/17S_F_E6_.fasta', 'merge_fixed/18D_F_C3_.fasta', 'merge_fixed/18S_F_F6_.fasta', 'merge_fixed/19D_G_B2_.fasta', 'merge_fixed/19S_G_H5_.fasta', 'merge_fixed/1D_F_A1_.fasta', 'merge_fixed/1S_F_A4_.fasta', 'merge_fixed/20D_F_D3_.fasta', 'merge_fixed/20S_F_C6_.fasta', 'merge_fixed/21D_F_E3_.fasta', 'merge_fixed/21S_F_H6_.fasta', 'merge_fixed/22D_G_G3_.fasta', 'merge_fixed/22S_G_A6_.fasta', 'merge_fixed/23D_F_F3_.fasta', 'merge_fixed/23S_F_D6_.fasta', 'merge_fixed/24D_G_H3_.fasta', 'merge_fixed/24S_G_B6_.fasta', 'merge_fixed/2D_F_B1_.fasta', 'merge_fixed/2S_F_B4_.fasta', 'merge_fixed/3D_F_C1_.fasta', 'merge_fixed/3S_F_C4_.fasta', 'merge_fixed/4D_G_D1_.fasta', 'merge_fixed/4S_G_D4_.fasta', 'merge_fixed/5D_G_E1_.fasta', 'merge_fixed/5S_G_E4_.fasta', 'merge_fixed/6D_G_F1_.fasta', 'merge_fixed/6S_G_F4_.fasta', 'merge_fixed/7D_G_G1_.fasta', 'merge_fixed/7S_G_H4_.fasta', 'merge_fixed/8D_G_H1_.fasta', 'merge_fixed/8S_G_G4_.fasta', 'merge_fixed/9D_G_A2_.fasta', 'merge_fixed/9S_G_A5_.fasta', 'merge_fixed/N_DOM_REPS_A7_.fasta', 'merge_fixed/N_GRA_A7_.fasta', 
#        '-L', 'merge/T4.fastq', 'merge/T6.fastq', 'merge/T7.fastq', 'merge/T8.fastq', 'merge/T9.fastq', 'merge/T10.fastq', 'merge/T11.fastq', 'merge/T12.fastq', 'merge/T13.fastq', 'merge/T14.fastq', 'merge/T15.fastq', 'merge/T16.fastq',
#        '-S', 'specifications.txt', 
#        '-o', 'numtdumper/', 
#        '-t', '4',
#        '-l', '418', 
#        '-p', '0', 
#        '-s', '5', 
#        '--matchpercent', '99.5'#,
#    #    '-T', 'numtdumper/5_denoise_coleoptera_fftnsi_UPGMA.nwk'
#        ])
#    
    args = parser.parse_args()
    
    # Find the file name
    
    filename = os.path.splitext(os.path.basename(args.asvs))[0]
    
    # Make the output directory
    
    if not os.path.exists(args.outputdirectory):
        os.makedirs(args.outputdirectory)
    
    ###########################################
    # READ AND PARSE FILTERING SPECIFICATIONS #
    ###########################################
    
    specs, termorder, thresholdcombos = assessmentcore.parse_specfile(args)
    
    # Resolve length specifications into minimum and maximum allowed
    
    minmaxbp = filterlength.resolve_length_spec(args)
    
    ##################
    # REPORT TO USER #
    ##################
    
    sys.stdout.write("\nWelcome to NUMTdumper, let's dump those NUMTs!\n\n"
                     f"Parsed {len(specs)} specification term"
                     f"{'s' if len(specs) > 1 else ''}, "
                     f"{len(thresholdcombos)} total threshold combination"
                     f"{'s' if len(thresholdcombos) > 1 else ''}\n")
    
    ###############
    # FIND CLADES #
    ###############
    #TODO: error catch for duplicate headers?
    clades, raw = findclades.find_clades(args, filename)
    
    #############
    # READ TAXA #
    #############
    
    taxa = None
    
    if args.taxgroups:
        sys.stdout.write("Reading taxa data\n")
        taxa = findtaxa.parse_taxa(args.taxgroups, raw['asvs'].keys())
    else:
        taxa = findtaxa.dummy_taxa(raw['asvs'].keys())
    
    ########################################
    # COMPUTE LIBRARY AND TOTAL READCOUNTS #
    ########################################
    
    sys.stdout.write("Matching library reads to ASVs to generate library ASV "
                     "counts.\n")
    
    librarycounts, totalcounts = findlibraries.count_asvs_in_libraries(
                                                               raw['asvs'],
                                                               args.libraries)
    
    # Output csv of library counts
    categorycounting.write_count_dict(librarycounts, raw['asvs'].keys(),
                                      os.path.join(args.outputdirectory,
                                                   f"{filename}_ASVcounts.csv"
                                                   ))
    
    ##########################
    # DESIGNATE CONTROL SETS #
    ##########################
    
    target, nontarget = assessmentcore.get_validated(raw, minmaxbp, args,
                                                     filename)
    
    ####################
    # CONSOLIDATE DATA #
    ####################
    
    data = {"total" : totalcounts,
           "library": librarycounts,
           "clade"  : clades,
           "taxon"  : taxa}
    
    ###########################################
    # GENERATE COUNTS AND SCORES, THEN ASSESS #
    ###########################################
    
    # Generate category counts for each specification
    sys.stdout.write("Generating binned counts\n")
    
    counts = assessmentcore.counts_from_spec(specs, data)
    
    # Calculate score for threshold combination
    sys.stdout.write("Assessing counts and scoring for each threshold "
                     "combination.\n")
    
    chunksize = math.ceil(len(thresholdcombos)/args.threads)
    with Pool(processes = args.threads) as pool:
        stats = pool.map(partial(assessmentcore.calc_stats, counts, termorder,
                                 set(raw['asvs'].keys()), target, nontarget,
                                 args.anythreshold, "standardised"),
                         thresholdcombos, chunksize)
    
#    stats = []
#    for c in thresholdcombos:
#        stats.append(assessmentcore.calc_stats(counts,
#                                               set(raw['asvs'].keys()),
#                                               target, nontarget,
#                                               args.anythreshold,
#                                               "standardised", c))
#    
    
#    
#    
#    counts_bak = counts
#    # Calculate score for threshold combination
#    
#    sys.stdout.write("Assessing counts and scoring for each threshold "
#                     "combination.\n")
#    
#    stats_filename = "statsdumptemp.pydata"
#    
#    try:
#        with open(stats_filename, "rb") as file:
#            stats = pickle.load(file)
#    except:
#        chunksize = math.floor(len(threshold_combinations)/args.threads)
#        with Pool(processes = args.threads) as pool:
#            stats = pool.map(partial(assessmentcore.calc_stats, counts, set(rawasvs.keys()), target, nontarget, args.anythreshold, "standardised" ),
#                         threshold_combinations, chunksize)
#        with open(stats_filename, "wb") as f:
#            pickle.dump(stats, f, pickle.HIGHEST_PROTOCOL)
#    
    
    # Analyse scores to find best sets
    
    #TODO: collapse different thresholds with identical outputs.
    
    sys.stdout.write("Identifying optimal threshold sets\n")
    
    minscore, minthresh, outtext = assessmentcore.get_minimum_thresholds(
                                                              stats, termorder,
                                                              thresholdcombos,
                                                              specs)
    
    #TODO: what to do if very large numbers of outputs with minimum score? Suppress outtext and output filtered ASVs for only a subset?
    
    #print(outtext)
    
    ##################
    # OUTPUT RESULTS #
    ##################
    
    # Output filtered haplotypes for threshold combination(s)
    
    sys.stdout.write("Writing filtered ASVs\n")
    
    assessmentcore.output_filtered_haplotypes(counts, minthresh, termorder,
                                              args.anythreshold, target, 
                                              nontarget, raw['path'], filename,
                                              args.outputdirectory)
    
    # Output thresholds and scores
    
    sys.stdout.write("Writing all data\n")
    
    assessmentcore.write_specs_and_stats(termorder, thresholdcombos, stats,
                                         os.path.join(args.outputdirectory,
                             f"{filename}_thresholds_scores_asvcounts.csv"))
    
    sys.stdout.write("\nNUMTs: dumped\n\n")
