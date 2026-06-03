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

from collections import defaultdict
from functools import partial
from shutil import which
from Bio import SeqIO

from . import core
from . import binning
from . import filterlength
from . import __version__
import csv


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
    tools = ['mafft', 'Rscript', 'bbmap.sh']
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
    parser.add_argument("--version", action="version", version=f"metaMATE {__version__}")
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
    corearg.add_argument("-o", "--output",
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
    corearg.add_argument("--enforcevalidation",
                         help="automatically remove all validated non-authentic ASVs and retain "
                              "all validated authentic ASVs in the final output, regardless of "
                              "threshold-based filtering (default: enabled). Use "
                              "--no-enforcevalidation to disable and output only threshold-based "
                              "filtering results",
                         action=argparse.BooleanOptionalAction,
                         default=True)

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

    # OTU binning variables
    otubin = coreparser.add_argument_group('otu binning')
    otubin.add_argument("--uc",
                        help="path to a .uc file (vsearch output) mapping ASVs to OTUs",
                        metavar="path", type=str)
    otubin.add_argument("--otu_fasta",
                        help="path to a fasta file containing OTU centroid sequences",
                        metavar="path", type=str)
    otubin.add_argument("--otu_table",
                        help="path to a table (CSV/TSV) containing OTU counts per library",
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
    refmatch.add_argument("--references2",
                          help="path to an optional second fasta of reference sequences "
                               "(e.g. a smaller, more local reference that better matches the "
                               "query sequences); it is concatenated with -R/--references and "
                               "the combined set is used for matching",
                          metavar="path", type=str)
    refmatch.add_argument("--refmatchlength",
                          help="the minimum alignment length to consider a match when "
                               "comparing ASVs against reference sequences (default is 80%% of "
                               "[calculated value of] -n/--minimumlength)",
                          type=int, metavar="n")
    # refmatch.add_argument("--refmatchpercent",
    #                       help="the minimum percent identity to consider a match when "
    #                            "comparing ASVs against reference sequences"
    #                            "[ CURRENTLY ONLY EXACT MATCHES ARE SUPPORTED ]",
    #                       type=float, default=99.9,
    #                       action=Range, minimum=0, maximum=100)
    refmatch.add_argument("--ignoreambigASVs",
                          help="ASVs that match the same reference will not be considered vaASV "
                          "(default is to count the most abundant one as verified authentic).",
                          action="store_true", default=False)
    refmatch.add_argument("--keeptemporaryfiles",
                          help="don't delete the temporary bbmap result "
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

    # Set up adaptive filter subparser
    adaptparser = subparsers.add_parser("filter-adaptive", parents=[coreparser],
                                        help="filter per-sample based on distribution of non-authentic ASVs")
    adaptparser._optionals.title = "adaptive filtering arguments"
    
    # Validation args (same as findparser) - duplicating necessary groups/args
    # Reusing findparser groups directly would be cleaner if possible but they are bound to findparser
    # So we re-add the necessary arguments.
    
    # Reference matching parameters
    refmatch = adaptparser.add_argument_group("reference-matching-based target "
                                              "identification")
    refmatch.add_argument("-R", "--references",
                          help="path to a fasta of known correct reference sequences",
                          metavar="path", type=str)
    refmatch.add_argument("--references2",
                          help="path to an optional second fasta of reference sequences "
                               "(e.g. a smaller, more local reference that better matches the "
                               "query sequences); it is concatenated with -R/--references and "
                               "the combined set is used for matching",
                          metavar="path", type=str)
    refmatch.add_argument("--refmatchlength",
                          help="the minimum alignment length to consider a match when "
                               "comparing ASVs against reference sequences (default is 80%% of "
                               "[calculated value of] -n/--minimumlength)",
                          type=int, metavar="n")
    refmatch.add_argument("--ignoreambigASVs",
                          help="ASVs that match the same reference will not be considered vaASV "
                          "(default is to count the most abundant one as verified authentic).",
                          action="store_true", default=False)
    refmatch.add_argument("--keeptemporaryfiles",
                          help="don't delete the temporary bbmap result "
                          "files generated during reference matching",
                          action="store_true", default=False)

    # Length parameters
    lengths = adaptparser.add_argument_group("length-based non-target identification")
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
    transl = adaptparser.add_argument_group("translation-based non-target identification")
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

    # Adaptive specific args
    adaptparser.add_argument("--percentile",
                             help="percentile of non-authentic ASVs to filter out (default 0.95)",
                             type=float, default=0.95,
                             action=Range, minimum=0, maximum=1)
    adaptparser.add_argument("--criteria",
                             help="criteria for filtering: 'verified_removed' (default) or 'estimated_removed'",
                             choices=['verified_removed', 'estimated_removed'],
                             default='verified_removed')

    args = parser.parse_args(arglist) if arglist else parser.parse_args()


    # Check for all required variables
    if args.mode in ['find', 'filter-adaptive']:

        # Ensure a value is supplied to libraries
        otu_args = [args.uc, args.otu_fasta, args.otu_table]
        otu_mode = any(otu_args)
        args.otu_mode = otu_mode
        if otu_mode and not all(otu_args):
            parser.error("if any of --uc, --otu_fasta, or --otu_table are provided, ALL must be provided")
            
        if not otu_mode and ((not args.libraries and not args.readmap)
                or (args.libraries and args.readmap)):
            parser.error("one and only one of -L/--libraries or -M/--readmap is required for "
                         "error finding (unless running in OTU mode)")
        
        # Only check metric if in find mode
        if args.mode == 'find' and args.scoremetric not in ['accuracy', 'precision', 'recall']:
            parser.error("-q/--scoremetric must be one of 'accuracy', 'precision' or 'recall'")
            
        # Ensure at least one reference is supplied
        if not args.references:
            parser.error("-R/--references is required for "
                         "error finding")
        # Check the length specification
        args, lset = filterlength.resolve_length_spec(args, parser)
        if not lset:
            parser.error("supply some length-based non-target identification specifications")
        # Set the matchlength defaults
        if not args.refmatchlength and args.references:
            args.refmatchlength = int(0.8 * args.minimumlength)
        # create output folder if it does not exist yet
        if not os.path.exists(args.output):
            os.makedirs(args.output)
 
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

        if args.resultcache:
            cache_otu_mode = core.check_cache_otu_mode(args.resultcache)
            otu_args = [args.uc, args.otu_fasta, args.otu_table]
            
            if cache_otu_mode:
                if not all(otu_args):
                     parser.error("The provided result cache was created in OTU mode. You MUST provide the OTU arguments (`--uc`, `--otu_fasta`, `--otu_table`) to dump these results correctly.")
            else:
                 # Warn if OTU args provided for non-OTU cache?
                 pass

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
            base_path = args.output
            if not os.path.splitext(base_path)[1]:
                 base_path = os.path.join(args.output, "dump")
            paths = core.make_resultset_paths(base_path, args.resultindex).values()
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


def subset_table(fasta_path, table_path, output_path, filtered_library_counts=None):
    """Filter an abundance table to retain only rows whose
    first-column ID appears in the given FASTA file.

    If `filtered_library_counts` is provided (mapping sample -> {id: count}),
    additionally zero out any cell whose value did not survive per-sample
    threshold filtering. This assumes the table has IDs as rows and samples
    as columns (header[0] is the row-ID column, header[1:] are sample names).
    """
    # Parse FASTA sequence IDs
    fasta_ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                fasta_ids.add(line.strip().lstrip(">"))

    sys.stdout.write(f"Subsetting table: {len(fasta_ids)} sequences in FASTA\n")

    # Detect separator
    with open(table_path) as f:
        firstline = f.readline().strip()
    if '\t' in firstline:
        sep = '\t'
    elif ',' in firstline:
        sep = ','
    else:
        sys.stderr.write(f"Warning: cannot detect separator in {table_path}, assuming tab\n")
        sep = '\t'

    # Read table, filter rows by first column
    with open(table_path, newline="") as f:
        reader = csv.reader(f, delimiter=sep)
        header = next(reader)
        kept_rows = []
        skipped = 0
        for row in reader:
            if row[0] in fasta_ids:
                kept_rows.append(row)
            else:
                skipped += 1

    sys.stdout.write(f"  OTUs kept: {len(kept_rows)}, removed: {skipped}\n")

    # Optionally zero out cells eliminated by per-sample threshold filter.
    cells_zeroed = 0
    if filtered_library_counts is not None:
        sample_names = header[1:]
        missing_samples = [s for s in sample_names if s not in filtered_library_counts]
        if missing_samples:
            sys.stderr.write(
                f"  Warning: header sample names not found in filtered counts "
                f"({', '.join(missing_samples)}); their columns will be zeroed.\n"
            )
        for row in kept_rows:
            row_id = row[0]
            for j, sample in enumerate(sample_names, start=1):
                kept_count = filtered_library_counts.get(sample, {}).get(row_id, 0)
                if kept_count == 0 and row[j] not in ("", "0"):
                    row[j] = "0"
                    cells_zeroed += 1
        sys.stdout.write(f"  cells zeroed by per-sample threshold: {cells_zeroed}\n")

    # Write filtered table using same separator
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter=sep)
        writer.writerow(header)
        for row in kept_rows:
            writer.writerow(row)

    sys.stdout.write(f"  Filtered table written to: {output_path}\n")


def main():
    #######################
    # INITIAL PREPARATION #
    #######################

    # Check for required programs
    check_tools()

    # Get command line arguments
    args = getcliargs()

    # Check for bad arguments
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
        
        # Check for OTU mode
        otu_args = [args.uc, args.otu_fasta, args.otu_table]
        otu_mode = any(otu_args)
        
        if otu_mode:
            if not all(otu_args):
                print("Error: In dump mode, if any OTU arguments are provided (`--uc`, `--otu_fasta`, `--otu_table`), ALL must be provided to output OTU sequences.")
                sys.exit(1)
            sys.stdout.write("Dump mode: OTU filtering enabled. Loading OTU sequences.\n")
            raw = {}
            # Load OTU sequences directly
            raw['asvs'] = SeqIO.to_dict(SeqIO.parse(args.otu_fasta, "fasta"))
            raw['path'] = args.otu_fasta
            aligned = {} # Not used in dump
        else:
            raw, aligned = binning.parse_asvs(args, True, '', tempbasepath)

        base_path = args.output
        if not os.path.splitext(base_path)[1]:
             base_path = os.path.join(args.output, "dump")
        core.write_resultset_asvs(set(raw['asvs'].keys()), args.resultindex, args.resultcache,
                                  raw['path'], base_path, args.mode)

        # Post-process: enforce validation on the output FASTAs
        if args.enforcevalidation:
            cache_dir = os.path.dirname(os.path.abspath(args.resultcache))
            asv_basename = os.path.splitext(os.path.basename(args.asvs))[0]
            paths = core.make_resultset_paths(base_path, args.resultindex)

            if otu_mode:
                otu_summary_path = os.path.join(cache_dir, "otu_summary.csv")
                if os.path.exists(otu_summary_path):
                    ev_target, ev_nontarget = core.load_validation_from_otu_summary(otu_summary_path)
                    for fasta_path in paths.values():
                        if os.path.exists(fasta_path):
                            core.enforce_validation_fasta(fasta_path, raw['path'],
                                                          ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic OTUs "
                                     f"included, {len(ev_nontarget)} non-authentic OTUs removed\n")
                else:
                    sys.stdout.write(f"Warning: OTU summary not found at {otu_summary_path}, "
                                     f"skipping validation enforcement\n")
            else:
                controlpath = os.path.join(cache_dir, asv_basename + "_control.txt")
                if os.path.exists(controlpath):
                    ev_target, ev_nontarget = core.load_validation_from_control(controlpath)
                    for fasta_path in paths.values():
                        if os.path.exists(fasta_path):
                            core.enforce_validation_fasta(fasta_path, raw['path'],
                                                          ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic ASVs "
                                     f"included, {len(ev_nontarget)} non-authentic ASVs removed\n")
                else:
                    sys.stdout.write(f"Warning: control file not found at {controlpath}, "
                                     f"skipping validation enforcement\n")

            # Subset the original table to match the filtered FASTA(s).
            # In OTU mode the FASTAs contain OTU centroid IDs, so subset the OTU
            # table; otherwise subset the ASV-level readmap.
            original_table = args.otu_table if otu_mode else args.readmap
            if original_table:
                for fasta_path in paths.values():
                    if os.path.exists(fasta_path):
                        table_out = os.path.splitext(fasta_path)[0] + "_table_filtered.txt"
                        subset_table(fasta_path, original_table, table_out)

        if os.path.exists(tempbasepath + "unaligned.fasta"):
            os.remove(tempbasepath + "unaligned.fasta")
        sys.stdout.write("\nCompleted dump\n\n")
        sys.stdout.flush()
        sys.exit(0)

    ###########################################
    # READ AND PARSE FILTERING SPECIFICATIONS #
    ###########################################

    specs, terms, nterm, nthresh, thresholds = core.parse_specs(args, 0)
    specs['otu_mode'] = args.otu_mode
    if args.mode != 'filter-adaptive':
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
        if args.otu_mode:
             sys.stdout.write("OTU mode active!\n")
             counts = ({}, {'total': {asv: 1 for asv in raw['asvs']}})
        else:
             sys.exit("Error: no library read count information supplied")

    librarycounts, totalcounts = counts
    librarysizes = [len(asvs) for lib, asvs in librarycounts.items()]
    librarysizes = {str(n): sum([1 for s in librarysizes if s == n]) for n in set(librarysizes)}
    try:
        pc1 = librarysizes['1'] / sum(librarysizes.values())
    except KeyError:
        pc1 = 0
    if pc1 >= 0.5:
        sys.stderr.write(f"Warning: {round(pc1 * 100, 0)}% of libraries only have one ASV!\n")

    # Output csv of library counts
    if args.libraries or args.readmap:
        core.write_count_dict(librarycounts, raw['asvs'].keys(), libcountpath)

    ##########################
    # DESIGNATE CONTROL SETS #
    ##########################

    target, nontarget = [set(), set()]
    if args.mode in ['find', 'filter-adaptive']:
        target, nontarget = core.get_validated(raw, args, baseinpath, totalcounts)

    ####################
    # APPLY OTU FILTER #
    ####################

    if args.otu_mode:
        sys.stdout.write("Applying OTU-level filtering using external OTU data...\n")
        
        # 1. Parse .uc file
        sys.stdout.write(f"Parsing UC file: {args.uc}\n")
        otus = binning.parse_uc(args.uc)
        
        # 2. Classify OTUs
        otu_to_asvs = defaultdict(list)
        for asv, otu in otus.items():
            otu_to_asvs[otu].append(asv)
            
        new_target = set()
        new_nontarget = set()
        asv_summary = []
        
        for otu, asv_list in otu_to_asvs.items():
            is_authentic = any(asv in target for asv in asv_list)
            is_nonauthentic = all(asv in nontarget for asv in asv_list)
            
            otu_status = "Unclassified"
            if is_authentic:
                new_target.add(otu)
                otu_status = "Authentic"
            elif is_nonauthentic:
                new_nontarget.add(otu)
                otu_status = "Non-Authentic"
            
            for asv in asv_list:
                asv_status = "Unclassified"
                if asv in target:
                    asv_status = "RefMatch"
                elif asv in nontarget:
                    asv_status = "NonAuthentic"
                asv_summary.append([asv, otu, asv_status, otu_status])
                
        # 3. Switch Data
        sys.stdout.write(f"Loading OTU sequences from {args.otu_fasta}\n")
        raw['asvs'] = SeqIO.to_dict(SeqIO.parse(args.otu_fasta, "fasta"))
        raw['path'] = args.otu_fasta
        
        sys.stdout.write(f"Loading OTU counts from {args.otu_table}\n")
        librarycounts, totalcounts = binning.parse_readmap(raw['asvs'], args.otu_table)
        
        # Update target/nontarget
        target = new_target
        nontarget = new_nontarget
        
        # Reset clades and taxa for OTUs (as they were calculated for ASVs)
        clades = binning.dummy_grouping(raw['asvs'].keys())
        taxa = binning.dummy_grouping(raw['asvs'].keys())
        
        # Write Summary File
        summary_path = os.path.join(args.output, "otu_summary.csv")
        with open(summary_path, 'w') as f:
            f.write("ASV,OTU,ASV_Status,OTU_Status\n")
            for row in asv_summary:
                f.write(",".join(row) + "\n")
        sys.stdout.write(f"Written OTU summary to {summary_path}\n")
        
        sys.stdout.write(f"OTU Filtering Complete. Proceeding with {len(raw['asvs'])} OTUs.\n")
        sys.stdout.write(f"  Authentic OTUs: {len(target)}\n")
        sys.stdout.write(f"  Non-Authentic OTUs: {len(nontarget)}\n")
        
        


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
    counts = {}
    if args.mode != 'filter-adaptive':
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

        dump_fasta = os.path.join(args.output, "dump_filtered.fasta")
        core.write_retained_asvs(raw['path'], dump_fasta, rejects)

        # Post-process: enforce validation if control/summary file is available
        if args.enforcevalidation:
            if hasattr(args, 'otu_mode') and args.otu_mode:
                otu_summary_path = os.path.join(args.output, "otu_summary.csv")
                if os.path.exists(otu_summary_path):
                    ev_target, ev_nontarget = core.load_validation_from_otu_summary(otu_summary_path)
                    core.enforce_validation_fasta(dump_fasta, raw['path'],
                                                  ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic OTUs "
                                     f"included, {len(ev_nontarget)} non-authentic OTUs removed\n")
                else:
                    sys.stdout.write(f"Warning: OTU summary not found at {otu_summary_path}, "
                                     f"skipping validation enforcement\n")
            else:
                controlpath = baseinpath + "_control.txt"
                if os.path.exists(controlpath):
                    ev_target, ev_nontarget = core.load_validation_from_control(controlpath)
                    core.enforce_validation_fasta(dump_fasta, raw['path'],
                                                  ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic ASVs "
                                     f"included, {len(ev_nontarget)} non-authentic ASVs removed\n")
                else:
                    sys.stdout.write(f"Warning: control file not found at {controlpath}, "
                                     f"skipping validation enforcement\n")

            # Subset the original table to match the filtered FASTA. In OTU mode,
            # subset the OTU table; otherwise use the ASV-level readmap (see note in
            # the filter-adaptive branch).
            if args.otu_mode:
                original_table = args.otu_table
            else:
                original_table = args.readmap if args.readmap else libcountpath
            if os.path.exists(original_table):
                table_out = os.path.splitext(dump_fasta)[0] + "_table_filtered.txt"
                subset_table(dump_fasta, original_table, table_out)
            else:
                sys.stderr.write(f"Warning: count table not found at {original_table}, "
                                 f"skipping table subsetting.\n")

        sys.stdout.write("\nCompleted dump\n\n")

    elif args.mode == 'filter-adaptive':
        sys.stdout.write(f"Running adaptive filtering (percentile: {args.percentile}, criteria: {args.criteria})\n")
        
        # adaptive filtering doesn't use the standard specs parsing/counts logic
        # instead it operates on the library counts directly
        
        output_csv = os.path.join(args.output, "filter-adaptive_table.csv")
        output_fasta = os.path.join(args.output, "filter-adaptive.fasta")
        output_summary = os.path.join(args.output, "filter-adaptive_summary.csv")
        
        filtered_library_counts, sample_thresholds, fallback_threshold, retained_asvs = core.adaptive_filter(
            raw['asvs'], librarycounts, target, nontarget,
            args.percentile, args.criteria,
            output_fasta, output_summary,
            enforcevalidation=args.enforcevalidation)

        # Post-process: enforce validation on output FASTA
        if args.enforcevalidation:
            ev_nontarget = set()
            if args.otu_mode:
                otu_summary_path = os.path.join(args.output, "otu_summary.csv")
                if os.path.exists(otu_summary_path):
                    ev_target, ev_nontarget = core.load_validation_from_otu_summary(otu_summary_path)
                    core.enforce_validation_fasta(output_fasta, raw['path'],
                                                  ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic OTUs "
                                     f"included, {len(ev_nontarget)} non-authentic OTUs removed\n")
                else:
                    sys.stdout.write(f"Warning: OTU summary not found at {otu_summary_path}, "
                                     f"skipping validation enforcement\n")
            else:
                controlpath = baseinpath + "_control.txt"
                if os.path.exists(controlpath):
                    ev_target, ev_nontarget = core.load_validation_from_control(controlpath)
                    core.enforce_validation_fasta(output_fasta, raw['path'],
                                                  ev_target, ev_nontarget)
                    sys.stdout.write(f"Enforced validation: {len(ev_target)} authentic ASVs "
                                     f"included, {len(ev_nontarget)} non-authentic ASVs removed\n")
                else:
                    sys.stdout.write(f"Warning: control file not found at {controlpath}, "
                                     f"skipping validation enforcement\n")

            core.write_adaptive_csv(filtered_library_counts, retained_asvs - ev_nontarget,
                                    sample_thresholds, fallback_threshold, output_csv)
        else:
            core.write_adaptive_csv(filtered_library_counts, retained_asvs,
                                    sample_thresholds, fallback_threshold, output_csv)

            # Subset the original table to match the filtered FASTA.
            # In OTU mode, subset the OTU table (rows are OTU centroid IDs); otherwise
            # use the ASV-level readmap. Falling back to args.readmap in OTU mode would
            # subset an ASV table by OTU centroid IDs, dropping all non-centroid member
            # ASV reads from each surviving OTU.
            if args.otu_mode:
                original_table = args.otu_table
            else:
                original_table = args.readmap if args.readmap else libcountpath
            if os.path.exists(original_table):
                table_out = os.path.splitext(output_fasta)[0] + "_table_filtered.txt"
                subset_table(output_fasta, original_table, table_out,
                             filtered_library_counts=filtered_library_counts)
            else:
                sys.stderr.write(f"Warning: count table not found at {original_table}, "
                                 f"skipping table subsetting.\n")

        sys.stdout.write(f"\nCompleted adaptive filter\n")
        sys.stdout.write(f"Abundance table written to: {output_csv}\n")
        sys.stdout.write(f"FASTA written to: {output_fasta}\n")
        sys.stdout.write(f"Summary table written to: {output_summary}\n\n")


# TODO add resume points?
if __name__ == "__main__":
    main()
    exit()
