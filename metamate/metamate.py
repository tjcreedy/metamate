#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" """

# Imports

import argparse
import os
import sys
import math
import textwrap as _textwrap
import multiprocessing

from functools import partial
from shutil import which

from metamate import core
from metamate import binning
from metamate import filterlength


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

def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def check_tools():
    """Check to see whether the required external tools are available."""
    tools = ['mafft', 'blastn', 'Rscript']
    for tool in tools:
        # tool = tools[1]
        if not is_tool(tool):
            err = (f"Error: cannot locate executable \'{tool}\'. Please ensure it is installed "
                   f"and available for command line execution.")
            sys.exit(err)


def getcliargs(arglist=None):
    # B D E F H I J K M N O P Q U V W Y Z
    # e f g h i j k m q u v w z

    parser = argparse.ArgumentParser(
        description="arguments can be passed on the command line or in a file, one per line, and "
                    "the file specified as @args.txt on the commandline",
        fromfile_prefix_chars='@')
    parser._positionals.title = "mode to run"
    parser._optionals.title = "optional arguments"

    subparsers = parser.add_subparsers(help='', dest='mode')
    subparsers.required = True
    # Core inputs on a separate parser which acts as parent to subs

    coreparser = argparse.ArgumentParser(add_help=False)
    corearg = coreparser.add_argument_group('universal inputs and options')
    corearg.add_argument("-A", "--asvs",
                         help="path to a fasta of unique sequences to filter",
                         required=True, metavar="path", type=str)
    corearg.add_argument("-L", "--libraries",
                         help="path to fastx file(s) of individual libraries/discrete samples "
                              "from which ASVs were drawn, or a single fastx with "
                              "\';samplename=NAME;\' or \';barcodelabel=NAME;\' annotations in "
                              "headers.",
                         metavar="path", type=str, nargs='*')
    corearg.add_argument("-M", "--readmap",
                         help="path to a comma- or tab- separated tabular text file containing "
                              "read counts for ASVs by library",
                         metavar="path", type=str)
    corearg.add_argument("-y", "--anyfail",
                         help="reject ASVs when any incidences fail to meet a threshold (default "
                              "is all incidences)",
                         action="store_true", default=False)
    corearg.add_argument("-o", "--outputname",
                         help="the base directory/file name path to which intermediate and final "
                              "output data should be written, file extensions will be added as "
                              "necessary",
                         required=True, metavar="name")
    corearg.add_argument("--realign",
                         help="force (re)alignment of the input ASVs",
                         action="store_true", default=False)
    corearg.add_argument("--overwrite",
                         help="force overwriting of intermediate and/or final output(s)",
                         action="store_true", default=False)
    corearg.add_argument("-t", "--threads",
                         help="number of threads to use (default 1)",
                         default=1, metavar="n", type=int)

    # Clade finder variables
    cladbin = coreparser.add_argument_group('clade binning')
    cladbin.add_argument("-d", "--divergence",
                         help="divergence level to use for assigning clades (default 0.2)",
                         default=0.2, type=float,
                         action=Range, minimum=0, maximum=1)
    cladbin.add_argument("-T", "--tree",
                         help="path to an ultrametric tree of the ASVs",
                         metavar="path", type=str)
    cladbin.add_argument("--distancemodel",
                         help="substitution model for UPGMA tree estimation (passed to R "
                              "dist.dna, default F84)",
                         default="F84", type=str, metavar="MOD")
    taxbin = coreparser.add_argument_group('taxon binning')
    taxbin.add_argument("-G", "--taxgroups",
                        help="path to a two-column csv file specifying the taxon for each ASV",
                        metavar="path", type=str)

    # Set up finding subparser
    findparser = subparsers.add_parser("find", parents=[coreparser],
                                       help="find errors by supplying threshold ranges and "
                                            "control specifications")
    findparser._optionals.title = "find arguments"
    findparser.add_argument("-S", "--specification",
                            help="path to a text file detailing the read count binning strategy "
                                 "and thresholds. An example can be found in the github",
                            required=True, metavar="path", type=str)
    findparser.add_argument("-q", "--scoremetric",
                            help="which of the three calculated scoring metrics ('accuracy', "
                                 "'precision' or 'recall') to use for selecting the top scoring "
                                 "set(s) (default accuracy)",
                            default='accuracy', metavar='metric',
                            type=str)
    findparser.add_argument("-g", "--generateASVresults",
                            help="generate fasta files of retained ASVs for for threshold sets: "
                                 "if no value is given, generate for top scoring set(s) only; "
                                 "otherwise generate for the given proportion of top scoring sets "
                                 "(default 0: no ASVs output)",
                            default=0, const=True, nargs='?',
                            action=Range, minimum=0, maximum=1,
                            type=float)

    # Reference matching parameters
    refmatch = findparser.add_argument_group("reference-matching-based target "
                                             "identification")
    refmatch.add_argument("-R", "--references",
                          help="path to a fasta of known correct reference sequences",
                          metavar="path", type=str)
    refmatch.add_argument("-D", "--blastdb",
                          help="path to a blast database, such as nt",
                          metavar="path", type=str)
    refmatch.add_argument("--refmatchlength",
                          help="the minimum alignment length to consider a BLAST match when "
                               "comparing ASVs against reference sequences (default is 80%% of "
                               "[calculated value of] -n/--minimumlength)",
                          type=int, metavar="n")
    refmatch.add_argument("--refmatchpercent",
                          help="the minimum percent identity to consider a BLAST match when "
                               "comparing ASVs against reference sequences",
                          type=float, default=99.9,
                          action=Range, minimum=0, maximum=100)
    refmatch.add_argument("--dbmatchlength",
                          help="the minimum alignment length to consider a BLAST match when "
                               "comparing ASVs against a blast database (default is 80%% of "
                               "[calculated value of] -n/--minimumlength)",
                          type=int, metavar="n")
    refmatch.add_argument("--dbmatchpercent",
                          help="the minimum percent identity to consider a BLAST match when "
                               "comparing ASVs against a blast database",
                          type=float, default=100,
                          action=Range, minimum=0, maximum=100)
    refmatch.add_argument("--keeptemporaryfiles",
                          help="don't delete the temporary blast database and/or blast result "
                               "files generated during reference matching",
                          action="store_true", default=False)

    # Length parameters
    lengths = findparser.add_argument_group("length-based non-target identification")
    lengths.add_argument("-n", "--minimumlength",
                         help="designate ASVs that are shorter than this value as non-target",
                         type=int, default=0, metavar="n")
    lengths.add_argument("-x", "--maximumlength",
                         help="designate ASVs that are longer than this value as non-target",
                         type=int, default=float('Inf'), metavar="n")
    lengths.add_argument("-l", "--expectedlength",
                         help="the expected length of the sequences",
                         type=int, metavar="n")
    lengths.add_argument("-b", "--basesvariation",
                         help="the number of bases of variation from the expected length outside "
                              "which ASVs should be designated as non-target",
                         type=int, metavar="n")
    lengths.add_argument("-p", "--percentvariation",
                         help="the percentage variation from the expected length outside which "
                              "ASVs should be designated as non-target",
                         type=float,
                         action=Range, minimum=0, maximum=100)
    lengths.add_argument("-c", "--codonsvariation",
                         help="the number of codons of variation from the expected length outside "
                              "which ASVs should be designated as non-target",
                         type=int, metavar="n")
    lengths.add_argument("--onlyvarybycodon",
                         help="designate ASVs that fall within other length thresholds but do not "
                              "vary by a multiple of 3 bases from the expected length as "
                              "non-target",
                         action="store_true", default=False)

    # Translation parameters
    transl = findparser.add_argument_group("translation-based non-target identification")
    transl.add_argument("-s", "--table",
                        help="the number referring to the translation table to use for "
                             "translation filtering",
                        metavar="path", default=5)
    transl.add_argument("-r", "--readingframe",
                        help="coding frame of sequences, if known",
                        type=int, choices={1, 2, 3})
    transl.add_argument("--detectionconfidence",
                        help="confidence level for detection of reading frame (default 0.95)",
                        type=float, default=0.95,
                        action=Range, minimum=0, maximum=1)
    transl.add_argument("--detectionminstops",
                        help="minimum number of stops to encounter for detection (default 100, "
                             "may need to decrease for few input ASVs)",
                        type=int, default=100, metavar="n")

    # Set up dumping sub
    dumpparser = subparsers.add_parser("dump", parents=[coreparser],
                                       help="dump NUMTs found in a previous run or with fixed "
                                            "thresholds")
    dumpparser._optionals.title = "arguments"
    # TODO: make sure _resultcache is correct name
    dumpparser.add_argument("-C", "--resultcache",
                            help="path to the _resultcache file from a previous run",
                            metavar='path', type=str)
    dumpparser.add_argument("-i", "--resultindex",
                            help="one or more indices of result(s) from a previous run from which "
                                 "to generate filtered ASVs",
                            type=int, metavar='n', nargs='*')
    dumpparser.add_argument("-S", "--specification",
                            help="one or more [category(/ies); metric; threshold] strings "
                                 "denoting the (multiplicative) specification for dumping errors. "
                                 "If provided via CLI, each [] string must be quoted",
                            metavar="'[C(s);M;T]'", type=str, nargs='*')

    args = parser.parse_args(arglist) if arglist else parser.parse_args()


    # Check for all required variables
    if args.mode == 'find':

        # Ensure a value is supplied to libraries
        if ((not args.libraries and not args.readmap)
                or (args.libraries and args.readmap)):
            parser.error("one and only one of -L/--libraries or -M/--readmap is required for "
                         "error finding")
        # Ensure the metric supplied is valid
        if args.scoremetric not in ['accuracy', 'precision', 'recall']:
            parser.error("-q/--scoremetric must be one of 'accuracy', 'precision' or 'recall'")
        # Ensure at least one reference is supplied
        if not args.references and not args.blastdb:
            parser.error("at least one of -R/--references and/or -D/--blastdb is required for "
                         "error finding")
        # Check the length specification
        args, lset = filterlength.resolve_length_spec(args, parser)
        if not lset:
            parser.error("supply some length-based non-target identification specifications")
        # Set the matchlength defaults
        if not args.dbmatchlength and args.blastdb:
            args.dbmatchlength = int(0.8 * args.minimumlength)
        if not args.refmatchlength and args.references:
            args.refmatchlength = int(0.8 * args.minimumlength)

    elif args.mode == 'dump':
        ressum = sum([args.resultcache is not None,
                      args.resultindex is not None])
        if args.specification:
            if ressum > 0:
                parser.error("-S/--specification is not compatible with any of -C/--resultcache, "
                             "-i/--resultindex")
            if (not args.libraries and not args.readmap) or (args.libraries and args.readmap):
                parser.error("one and only one of -L/--libraries or -M/--readmap is required if "
                             "providing -S/--specification")
        elif ressum == 0:
            parser.error("either -S/--specification, or both of -C/--resultcache, "
                         "-i/--resultindex is required for dumping errors")
        elif ressum == 1:
            parser.error("both -C/--resultcache and -i/--resultindex are required if either is "
                         "specified")

        # Check outputs
        if args.mode == 'find' or args.specification:
            if os.path.exists(args.output):
                if not args.overwrite:
                    sys.stderr.write(f"The directory {args.output} for intermediate files exists, "
                                     f"will attempt to resume from any files within, to overwrite "
                                     f"instead set --overwrite\n")
            else:
                os.makedirs(args.output)
        else:  # Must be args.mode == 'dump' and --resultindex of 1 or more values
            paths = core.make_resultset_paths(args.output, args.resultindex).values()
            exists = [os.path.exists(p) for p in paths]
            if any(exists) and not args.overwrite:
                if len(args.resultindex) > 1:
                    parser.error(f"{sum(exists)}/{len(args.resultindex)} of the expected output "
                                 f"fastas for these result indices already exists, set --overwrite "
                                 f"to overwrite")
                else:
                    parser.error(f"the output fasta {paths[0]} already exists, set --overwrite to "
                                 f"overwrite")

    sys.stderr.flush()
    return args


def main():
    #######################
    # INITIAL PREPARATION #
    #######################

    # Check for required programs
    check_tools()

    # Get inputs
    args = getcliargs()
    # args = getcliargs('find '
    #                   '-A /home/thomas/QMRmeta/2_aln.fasta '
    #                   '-R /home/thomas/QMRmeta/907_Qinling_barcode_by_NAPselect.fasta '
    #                   '-T /home/thomas/QMRmeta/2_aln_UPGMA.nwk '
    #                   '-M /home/thomas/QMRmeta/2_asvtable.tsv '
    #                   '-S /home/thomas/QMRmeta/specifications.txt '
    #                   '-o /home/thomas/QMRmeta/3_metamate/ '
    #                   '-s 5 -l 418 -p 0 -t 20 '
    #                   '-D ~/db/NCBI/NT/nt '
    #                   '-G /home/thomas/QMRmeta/1_taxonomy.csv'.split(' '))
    # args = getcliargs('dump '
    #                   '-o /home/thomc/QMRmeta/3_metamate.fasta '
    #                   '-C /home/thomc/QMRmeta/3_metamate/2_aln_resultcache '
    #                   '-i 1011 '
    #                   '-A /home/thomc/QMRmeta/1_taxfilter.fasta'.split(' '))

    sys.stdout.write(f"\nWelcome to metaMATE, let's {args.mode}!\n\n")

    baseinpath = os.path.join(args.output, os.path.splitext(os.path.basename(args.asvs))[0])

    #################################################
    # DO EXTRACTION IF EXTRACTING FROM PREVIOUS RUN #
    #################################################

    if args.mode == 'dump' and args.specification is None:
        tempbasepath = os.path.join('.', 'asvtemp')
        raw, aligned = binning.parse_asvs(args, True, '', tempbasepath)
        core.write_resultset_asvs(set(raw['asvs'].keys()), args.resultindex, args.resultcache,
                                  raw['path'], args.output, args.mode)
        if os.path.exists(tempbasepath + "unaligned.fasta"):
            os.remove(tempbasepath + "unaligned.fasta")
        sys.stdout.write("\nCompleted dump\n\n")
        exit()

    ###########################################
    # READ AND PARSE FILTERING SPECIFICATIONS #
    ###########################################

    specs, terms, nterm, nthresh, thresholds = core.parse_specs(args, 0)

    sys.stdout.write(f"Parsed {nterm} additive specification term"
                     f"{'s' if nterm > 1 else ''}, comprising "
                     f"{len(specs['name'])} bin strateg"
                     f"{'ies' if len(specs['name']) > 1 else 'y'}")
    if args.mode == 'find':
        sys.stdout.write(f" and totalling {nthresh} unique threshold "
                         f"set{'s' if nthresh > 1 else ''}.\n")
    else:
        sys.stdout.write('.\n')

    ###############
    # FIND CLADES #
    ###############
    # TODO: error catch for duplicate headers?

    if 'clade' in terms:
        clades, raw = binning.find_clades(args, baseinpath)
    else:
        raw, aligned = binning.parse_asvs(args, True, '', baseinpath)
        clades = binning.dummy_grouping(raw['asvs'].keys())

    #############
    # READ TAXA #
    #############

    if args.taxgroups:
        sys.stdout.write("Reading taxa data\n")
        taxa = binning.parse_taxa(args.taxgroups, raw['asvs'].keys())
    else:
        taxa = binning.dummy_grouping(raw['asvs'].keys())

    ########################################
    # COMPUTE LIBRARY AND TOTAL READCOUNTS #
    ########################################

    libcountpath = baseinpath + "_ASVcounts.csv"

    if os.path.exists(libcountpath) and not args.overwrite:
        sys.stdout.write(f"Resuming with read mapping table of library ASV counts found at "
                         f"{libcountpath}.\n")
        counts = binning.parse_readmap(raw['asvs'], libcountpath)
    elif args.readmap:
        sys.stdout.write("Parsing read mapping table to generate library ASV counts.\n")
        counts = binning.parse_readmap(raw['asvs'], args.readmap)
    elif args.libraries:
        sys.stdout.write("Matching library reads to ASVs to generate library ASV counts.\n")
        counts = binning.count_asvs_in_libraries(raw['asvs'], args.libraries)
    else:
        sys.exit("Error: no library read count information supplied")

    librarycounts, totalcounts = counts

    librarysizes = [len(asvs) for lib, asvs in librarycounts.items()]
    librarysizes = {str(n): sum([1 for s in librarysizes if s == n]) for n in set(librarysizes)}
    pc1 = librarysizes['1'] / sum(librarysizes.values())
    if pc1 >= 0.5:
        sys.stderr.write(f"Warning: {round(pc1 * 100, 0)}% of libraries only have one ASV!\n")

    # Output csv of library counts
    if args.libraries or args.readmap:
        core.write_count_dict(librarycounts, raw['asvs'].keys(), libcountpath)

    ##########################
    # DESIGNATE CONTROL SETS #
    ##########################

    target, nontarget = [{}, {}]
    if args.mode == 'find':
        target, nontarget = core.get_validated(raw, args, baseinpath)

    ####################
    # CONSOLIDATE DATA #
    ####################

    data = {"total": totalcounts,
            "library": librarycounts,
            "clade": clades,
            "taxon": taxa}

    ###################
    # GENERATE COUNTS #
    ###################

    # Generate category counts for each specification
    sys.stdout.write("Generating binned counts\n")

    counts = core.counts_from_spec(specs, data)

    ####################
    # SCORE AND OUTPUT #
    ####################

    if args.mode == 'find':

        # SCORING

        # Initialise the queue manager and pool
        manager = multiprocessing.Manager()
        pool = multiprocessing.Pool(args.threads + 2)

        # Start writer
        queues, watchers = core.start_writers(pool, manager, specs, args.output, nthresh)

        # Calculate score for threshold combination
        sys.stdout.write("Assessing counts and scoring for each threshold combination.\n\n")

        chunksize = math.ceil(nthresh / args.threads)
        scores = pool.map(partial(core.assess_numts,
                                  specs['name'], counts, args.anyfail, set(raw['asvs'].keys()),
                                  target, nontarget, queues),
                          thresholds, chunksize)
        # scores = []
        # for threshvals in thresholds:
        #     out = core.assess_numts(specs['name'], counts, args.anyfail, set(raw['asvs'].keys()),
        #                             target, nontarget, False, threshvals)
        #     scores.append(out)

        # Close queues and pool
        for q in queues:
            q.put(None)
        pool.close()
        pool.join()

        # Get any errors in writers
        for w in watchers:
            w.get()

        # Parse scores
        scoresort = core.find_best_score(scores, args.scoremetric)

        # OUTPUT

        hashcachepath = os.path.join(args.output, "hashcache")

        # Output ASVs if requested
        if args.generateASVresults > 0:
            args.generateASVresults = 1
            resultcachepath = os.path.join(args.output, "resultcache")

            resultsets = core.generate_resultsets(args.generateASVresults, scoresort, hashcachepath)
            core.write_resultset_asvs(set(raw['asvs'].keys()), resultsets, resultcachepath,
                                      raw['path'], args.output, args.mode)

        # Output thresholds and scores
        os.remove(hashcachepath)
        sys.stdout.write("\nCompleted find\n\n")

    elif args.mode == 'dump':

        # SCORE

        sys.stdout.write("Assessing counts against thresholds\n")
        rejects = core.apply_reject(specs['name'], next(thresholds)[2:], counts, args.anyfail)

        # OUTPUT

        core.write_retained_asvs(raw['path'], args.output + ".fasta", rejects)
        sys.stdout.write("\nCompleted dump\n\n")


# TODO add resume points?
if __name__ == "__main__":
    main()
    exit()
