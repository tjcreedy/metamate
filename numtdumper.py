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

# Function definitions

def getcliargs(arglist = None):
    # B D E F H I J K M N O P Q U V W Y Z
    # e f g h i j k m q u v w z

    parser = argparse.ArgumentParser(description = "arguments can be passed on the command line or in a file, one per line, and the file specified as @args.txt on the commandline",
                                     fromfile_prefix_chars = '@')
    parser._positionals.title = "action to perform"
    parser._optionals.title = "optional arguments"

    subparsers = parser.add_subparsers(help = '', dest = 'action')
    subparsers.required = True
        # Core inputs on a separate parser which acts as parent to subs

    coreparser = argparse.ArgumentParser(add_help = False)
    core = coreparser.add_argument_group('universal inputs and options')
    core.add_argument("-A", "--asvs",
                      help = "path to a fasta of unique sequences to filter",
                      required = True, metavar = "path", type = str)
    core.add_argument("-L", "--libraries",
                      help = "path to fastx file(s) of individual libraries"
                             "/discrete samples from which ASVs were drawn, "
                             "or a single fastx with ;samplename=.*; or "
                             ";barcodelabel=.*; annotations in headers.",
                      metavar = "path", type = str, nargs = '*')
    core.add_argument("-y", "--anyfail",
                      help = "reject ASVs when any incidences fail to meet "
                             "a threshold (default is all incidences)",
                      action = "store_true", default = False)
    core.add_argument("-o", "--outputdirectory",
                      help = "output directory (default is current directory)",
                      default = os.getcwd(), metavar = "path/")
    core.add_argument("--realign",
                      help = "force (re)alignment of the input ASVs",
                      action = "store_true", default = False)
    core.add_argument("--overwrite",
                      help = "force overwriting of output file or directory "
                             "if it already exists",
                      action = "store_true", default = False)

        # Clade finder variables
    cladbin = coreparser.add_argument_group('clade binning')
    cladbin.add_argument("-d", "--divergence",
                         help = "divergence level to use for assigning "
                                "clades (default is 0.2)",
                         default = 0.2, type = float,
                         action = Range, minimum = 0, maximum = 1)
    cladbin.add_argument("-T", "--tree",
                         help = "path to an ultrametric tree of the ASVs",
                         metavar = "path", type = str)
    cladbin.add_argument("--distancemodel",
                         help = "substitution model for UPGMA tree "
                                "estimation (passed to R dist.dna, default "
                                "F84)",
                         default = "F84", type = str, metavar = "MOD")
    taxbin = coreparser.add_argument_group('taxon binning')
    taxbin.add_argument("-G", "--taxgroups",
                        help = "path to a two-column csv file specifying the "
                               "taxon for each ASV",
                        metavar = "path", type = str)

        #Set up finding subparser

    findparser = subparsers.add_parser("find", parents = [coreparser],
                                       help = "find NUMTs by supplying "
                                              "threshold ranges and control "
                                              "specifications")
    findparser._optionals.title = "find arguments"
    findparser.add_argument("-S", "--specification",
                            help = "path to a text file detailing the read "
                                   "count binning strategy and thresholds",
                            metavar = "path", type = str)
    findparser.add_argument("-g", "--generateASVresults",
                            help = "generate fasta files of retained ASVs for "
                                   "for threshold sets: if no value is given, "
                                   "generate for top scoring set(s) only; "
                                   "otherwise generate for the given "
                                   "proportion of top scoring sets (default "
                                   "0: no ASVs output)",
                            default = 0, const = True, nargs = '?',
                            action = Range, minimum = 0, maximum = 1,
                            type = float)
    findparser.add_argument("-t", "--threads",
                            help = "number of threads to use (default 1)",
                            default = 1, metavar = "n", type = int)

        # Reference matching parameters
    refmatch = findparser.add_argument_group("reference-matching-based target "
                                             "identification")
    refmatch.add_argument("-R", "--references",
                          help = "path to a fasta of known correct reference "
                                 "sequences",
                          metavar = "path", type = str)
    refmatch.add_argument("-D", "--blastdb",
                          help = "path to a blast database, such as nt",
                          metavar = "path", type = str)
    refmatch.add_argument("--refmatchlength",
                          help = "the minimum alignment length to consider a "
                                 "BLAST match when comparing ASVs against "
                                 "reference sequences (default is 80% of "
                                 "[calculated value of] -n/--minimumlength)",
                          type = int, metavar = "n")
    refmatch.add_argument("--refmatchpercent",
                          help = "the minimum percent identity to consider a "
                                 "BLAST match when comparing ASVs against "
                                 "reference sequences",
                          type = float, default = 99.9,
                          action = Range, minimum = 0, maximum = 100)
    refmatch.add_argument("--dbmatchlength",
                          help = "the minimum alignment length to consider a "
                                 "BLAST match when comparing ASVs against "
                                 "a blast database (default is 80% of "
                                 "[calculated value of] -n/--minimumlength)",
                          type = int, metavar = "n")
    refmatch.add_argument("--dbmatchpercent",
                          help = "the minimum percent identity to consider a "
                                 "BLAST match when comparing ASVs against "
                                 "a blast database",
                          type = float, default = 100,
                          action = Range, minimum = 0, maximum = 100)

        # Length parameters
    lengths = findparser.add_argument_group("length-based non-target "
                                            "identification")
    lengths.add_argument("-n", "--minimumlength",
                         help = "designate ASVs that are shorter than this "
                                "value as non-target",
                         type = int, default = 0, metavar = "n")
    lengths.add_argument("-x", "--maximumlength",
                         help = "designate ASVs that are longer than this "
                                "value as non-target",
                         type = int, default = float('Inf'), metavar = "n")
    lengths.add_argument("-l", "--expectedlength",
                         help = "the expected length of the sequences",
                         type = int, metavar = "n")
    lengths.add_argument("-b", "--basesvariation",
                         help = "the number of bases of variation from the "
                                "expected length outside which ASVs should "
                                "be designated as non-target",
                        type = int, metavar = "n")
    lengths.add_argument("-p", "--percentvariation",
                         help = "the percentage variation from the expected "
                                "length outside which ASVs should be "
                                "designated as non-target",
                         type = float,
                         action = Range, minimum = 0, maximum = 100)
    lengths.add_argument("-c", "--codonsvariation",
                         help = "the number of codons of variation from the "
                                "expected length outside which ASVs should be "
                                "designated as non-target",
                         type = int, metavar = "n")
    lengths.add_argument("--onlyvarybycodon",
                         help = "designate ASVs that fall within other length "
                                "thresholds but do not vary by a multiple of "
                                "3  bases from the expected length as "
                                "non-target",
                         action = "store_true", default = False)

        # Translation parameters
    transl = findparser.add_argument_group("translation-based non-target "
                                           "identification")
    transl.add_argument("-s", "--table",
                        help = "the number referring to the translation table "
                               "to use for translation filtering",
                        metavar = "path", default = 5)
    transl.add_argument("-r", "--readingframe",
                        help = "coding frame of sequences, if known",
                        type = int, choices = {1,2,3})
    transl.add_argument("--detectionconfidence",
                        help = "confidence level for detection of reading "
                               "frame (default 0.95)",
                        type = float, default = 0.95,
                        action = Range, minimum = 0, maximum = 1)
    transl.add_argument("--detectionminstops",
                        help = "minimum number of stops to encounter for "
                               "detection (default 100, may need to decrease "
                               "for few input ASVs)",
                        type = int, default = 100, metavar = "n")

        #Set up dumping sub
    dumpparser = subparsers.add_parser("dump", parents = [coreparser],
                                       help = "dump NUMTs found in a previous "
                                              "run or with fixed thresholds")
    dumpparser._optionals.title = "arguments"
    dumpparser.add_argument("-f", "--outfasta",
                            help = "output file name, if only retreiving "
                                   "result ASVs of a single index ",
                            required = False, metavar = "path")
    #TODO: make sure _resultcache is correct name
    dumpparser.add_argument("-C", "--resultcache",
                            help = "path to the _resultcache file from a "
                                   "previous run ",
                            metavar = 'path', type = str)
    dumpparser.add_argument("-i", "--resultindex",
                            help = "one or more indices of result(s) from a "
                                   "previous run from which to generate "
                                   "filtered ASVs",
                            type = int, metavar = 'n', nargs = '*')
    dumpparser.add_argument("-S", "--specification",
                            help = "one or more [category(/ies); metric; "
                                   "threshold] strings denoting the "
                                   "(multiplicative) specification for "
                                   "dumping NUMTs. If provided via CLI, "
                                   "each [] string must be quoted",
                            metavar = "'[C(s);M;T]'", type = str, nargs = '*')
    
    args = parser.parse_args(arglist) if arglist else parser.parse_args()
    
    # Check for all required variables
    
    if args.action == 'find':
        # Ensure a value is supplied to libraries
        if not args.specification :
            args.specification = os.path.join(os.path.dirname(__file__),
                                               'specifications.txt')
            sys.stderr.write( "Warning: -S/--specifications not supplied, "
                              "falling back to the default specifications in "
                             f"{args.specification}\n")
        if not args.libraries:
            parser.error('-L/--libraries is required for NUMT finding')
        # Ensure at least one reference is supplied
        if not args.references and not args.blastdb:
            parser.error('at least one of -R/--references and/or '
                         '-D/--blastdb is required for NUMT finding')
        # Check the length specification
        args, lset = filterlength.resolve_length_spec(args, parser)
        if not lset:
            parser.error('supply some length-based non-target identification '
                         'specifications')
        # Set the matchlength defaults
        if not args.dbmatchlength and args.blastdb:
            args.dbmatchlength = int(0.8 * args.minimumlength)
        if not args.refmatchlength and args.references:
            args.refmatchlength = int(0.8 * args.minimumlength)
        args.outfasta = None
    
    elif args.action == 'dump':
        ressum = sum([args.resultcache is not None,
                      args.resultindex is not None])
        if args.specification:
            if ressum > 0:
                parser.error('-S/--specification is not compatible with any of'
                             '-C/--resultcache, -i/--resultindex')
            if not args.libraries:
                parser.error('-L/--libraries is required if providing '
                             '-S/--specification')
        elif ressum == 0:
            parser.error('either -S/--specification, or both of '
                         '-C/--resultcache, -i/--resultindex is required for '
                         'NUMT dumping')
        elif ressum == 1:
            parser.error('both -C/--resultcache and -i/--resultindex are '
                         'required if either is specified')
        elif args.outfasta and len(args.resultindex) > 1:
            parser.error('-f/--outfasta is redundant when passing >1 indices '
                         'to -i/--resultindex. Supply -o/--outputdirectory '
                         'instead')
            args.outfasta = None
    
    if not args.overwrite:
        for a, t in zip([args.outputdirectory, args.outfasta],
                        ['-o/--outputdirectory', '-f/--outfasta']):
            if a and os.path.exists(a):
                sys.exit(f"{t} {a} exists but --overwrite is not set.")
    
    sys.stderr.flush()
    return(args)

kwargs = {
'al' : ['find',
            '-A', '6_coleoptera_fftnsi.fasta',
            '-R', 'dummy_reference.fasta',
#        #'-L', 'merge_fixed/10D_F_C2_.fasta', 'merge_fixed/10S_F_B5_.fasta', 'merge_fixed/11D_F_D2_.fasta', 'merge_fixed/11S_F_C5_.fasta', 'merge_fixed/12D_G_E2_.fasta', 'merge_fixed/12S_G_D5_.fasta', 'merge_fixed/13D_G_F2_.fasta', 'merge/13S_G_G6_.fasta', 'merge_fixed/14D_G_G2_.fasta', 'merge_fixed/14S_G_E5_.fasta', 'merge_fixed/15D_F_H2_.fasta', 'merge_fixed/15S_F_G5_.fasta', 'merge_fixed/16D_F_A3_.fasta', 'merge_fixed/16S_F_F5_.fasta', 'merge_fixed/17D_F_B3_.fasta', 'merge_fixed/17S_F_E6_.fasta', 'merge_fixed/18D_F_C3_.fasta', 'merge_fixed/18S_F_F6_.fasta', 'merge_fixed/19D_G_B2_.fasta', 'merge_fixed/19S_G_H5_.fasta', 'merge_fixed/1D_F_A1_.fasta', 'merge_fixed/1S_F_A4_.fasta', 'merge_fixed/20D_F_D3_.fasta', 'merge_fixed/20S_F_C6_.fasta', 'merge_fixed/21D_F_E3_.fasta', 'merge_fixed/21S_F_H6_.fasta', 'merge_fixed/22D_G_G3_.fasta', 'merge_fixed/22S_G_A6_.fasta', 'merge_fixed/23D_F_F3_.fasta', 'merge_fixed/23S_F_D6_.fasta', 'merge_fixed/24D_G_H3_.fasta', 'merge_fixed/24S_G_B6_.fasta', 'merge_fixed/2D_F_B1_.fasta', 'merge_fixed/2S_F_B4_.fasta', 'merge_fixed/3D_F_C1_.fasta', 'merge_fixed/3S_F_C4_.fasta', 'merge_fixed/4D_G_D1_.fasta', 'merge_fixed/4S_G_D4_.fasta', 'merge_fixed/5D_G_E1_.fasta', 'merge_fixed/5S_G_E4_.fasta', 'merge_fixed/6D_G_F1_.fasta', 'merge_fixed/6S_G_F4_.fasta', 'merge_fixed/7D_G_G1_.fasta', 'merge_fixed/7S_G_H4_.fasta', 'merge_fixed/8D_G_H1_.fasta', 'merge_fixed/8S_G_G4_.fasta', 'merge_fixed/9D_G_A2_.fasta', 'merge_fixed/9S_G_A5_.fasta', 'merge_fixed/N_DOM_REPS_A7_.fasta', 'merge_fixed/N_GRA_A7_.fasta',
        #'-L', 'merge/T4.fastq', 'merge/T6.fastq', 'merge/T7.fastq', 'merge/T8.fastq', 'merge/T9.fastq', 'merge/T10.fastq', 'merge/T11.fastq', 'merge/T12.fastq', 'merge/T13.fastq', 'merge/T14.fastq', 'merge/T15.fastq', 'merge/T16.fastq',
        '-L', '1_concat.fastq',
        '-S', '/home/thomas/Documents/programming/bioinformatics/numtdumper/specifications.txt',
        '-o', 'livetest/',
        '-l', '418',
        '-p', '0',
        '-s', '5',
        '-t', '2',
        #'--realign',
        '--refmatchpercent', '100'#,
        #'-g'#,
#    #    '-T', 'numtdumper/5_denoise_coleoptera_fftnsi_UPGMA.nwk'
        ],
# 'al' : ['dump', 
#         '-A', '6_coleoptera.fasta', 
#         '-C', 'test1/6_coleoptera_resultcache',
#         '-i', '43',
#         '-f', 'test8.fa'],
# 'al' : ['dump',
#         '-A', '6_coleoptera.fasta',
#         '-S', '[library; n; 3]', '[library; p; 0.0025]',
#         '[library|clade; p; 0.04]',
#         '-f', 'test10.fa'],
'sd' : "/home/thomas/Documents/programming/bioinformatics/numtdumper/",
'wd' : "/home/thomas/seqtesting/NUMTdumper/amm/"
}

def main(**kwargs):
    # Set debug variables
    if kwargs.get('wd', None):
        os.chdir(kwargs.get('wd'))
    if kwargs.get('sd', None):
        scriptdir = kwargs.get('sd')
    
    #######################
    # INITIAL PREPARATION #
    #######################
    
    # Get inputs
    args = getcliargs(kwargs.get('al')) if kwargs.get('al') else getcliargs()
    
    # Find the file name
    infilename = os.path.splitext(os.path.basename(args.asvs))[0]
    outfilename = infilename
    if args.action == 'dump' and args.outfasta:
        outfilename = os.path.splitext(os.path.basename(args.outfasta))[0]
    
    # Make the output directory
    
    if args.outputdirectory and not os.path.exists(args.outputdirectory):
        os.makedirs(args.outputdirectory)
    
    sys.stdout.write(f"\nWelcome to NUMTdumper, let's {args.action} "
                      "those NUMTs!\n\n")
    
    #################################################
    # DO EXTRACTION IF EXTRACTING FROM PREVIOUS RUN #
    #################################################
    
    if args.action == 'dump' and args.specification is None:
        
        raw, aligned = findclades.parse_asvs(args, False, '',
                                             os.path.join('.', 'asvtemp'))
        
        stores = assessmentcore.parse_resultcache(args.resultcache,
                                                  set(raw['asvs'].keys()))
        
        assessmentcore.write_resultset_asvs(set(raw['asvs'].keys()),
                                            outfilename, raw['path'],
                                            os.getcwd(),
                                            args.resultindex, stores,
                                            args.action, args.outfasta)
        
        
        sys.stdout.write("\nNUMTs: dumped\n\n")
        exit()
        
    ###########################################
    # READ AND PARSE FILTERING SPECIFICATIONS #
    ###########################################
    
    specs, nterm, nthresh, terms, thresholds = assessmentcore.parse_specs(
                                                                      args, 0)
    
    sys.stdout.write(f"Parsed {nterm} additive specification term"
                     f"{'s' if nterm > 1 else ''}, comprising "
                     f"{len(specs['name'])} bin strateg"
                     f"{'ies' if len(specs['name']) > 1 else 'y'}")
    if args.action == 'find':
        sys.stdout.write(f"and totalling {nthresh} unique threshold "
                         f"set{'s' if nthresh > 1 else ''}\n")
    else:
        sys.stdout.write('\n')
    
    ###############
    # FIND CLADES #
    ###############
    #TODO: error catch for duplicate headers?
    clades, raw = findclades.find_clades(args, infilename)
    
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
                                      f"{infilename}_ASVcounts.csv"))
    
    ##########################
    # DESIGNATE CONTROL SETS #
    ##########################
    
    if args.action == 'find':
        target, nontarget = assessmentcore.get_validated(raw, args,
                                                         infilename)
    
    ####################
    # CONSOLIDATE DATA #
    ####################
    
    data = {"total" : totalcounts,
           "library": librarycounts,
           "clade"  : clades,
           "taxon"  : taxa}
    
    ##############################
    # GENERATE COUNTS AND SCORES #
    ##############################
    
    # Generate category counts for each specification
    sys.stdout.write("Generating binned counts\n")
    
    counts = assessmentcore.counts_from_spec(specs, data)
    
    if args.action == 'find':
        
        # Calculate score for threshold combination
        sys.stdout.write("Assessing counts and scoring for each threshold "
                         "combination.\n")
        
        chunksize = math.ceil(nthresh/args.threads)
        with Pool(processes = args.threads) as pool:
            stats = pool.map(partial(assessmentcore.assess_numts,
                                     specs['name'], counts, args.anyfail,
                                     set(raw['asvs'].keys()), target,
                                     nontarget, "standardised", 0.5),
                             thresholds, chunksize)
        
        scoresort = assessmentcore.find_best_score(stats, len(specs['name']))
        
    elif args.action == 'dump':
        rejects = assessmentcore.apply_reject(specs['name'], next(thresholds),
                                              counts, args.anyfail)
    
    ##################
    # OUTPUT RESULTS #
    ##################
    
    if args.action == 'find':
        
        # Output ASVs if requested
        if args.generateASVresults > 0:
            resultsets = assessmentcore.generate_resultsets(
                                                       args.generateASVresults,
                                                       stats, scoresort,
                                                       len(specs['name']))
            sys.stdout.write(f"Writing {len(resultsets)} filtered ASV files\n")
            assessmentcore.write_resultset_asvs(set(raw['asvs'].keys()),
                                                outfilename, raw['path'],
                                                args.outputdirectory,
                                                resultsets, stats,
                                                args.action, None)
        
        # Output thresholds and scores
        sys.stdout.write("Writing statistics and result cache\n")
        assessmentcore.write_stats_and_cache(specs, stats, terms,
                                             infilename, args.outputdirectory)
        sys.stdout.write("\nNUMTs: found?\n\n")
    
    elif args.action == 'dump':
        assessmentcore.write_retained_asvs(raw['path'], args.outfasta, 
                                           rejects)
        sys.stdout.write("\nNUMTs: dumped\n\n")

# TODO add resume points?
if __name__ == "__main__":
    main()
    exit()