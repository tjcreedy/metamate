#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core functions for designating numts and optimising parameters"""

# Imports

import os
import sys
import shutil
import re
import gzip
import itertools
import time
import datetime

import numpy

from collections import defaultdict
from functools import partial, reduce

from Bio.SeqIO.FastaIO import SimpleFastaParser

from . import filterlength
from . import filtertranslate
from . import filterreference


def count_categories(catdict, metric):
    # catdict, metric = [data[specdict['terms'][0]], specdict['metric']]

    # Set up dict of sets for scores
    counts = defaultdict(set)

    # Work through each category
    for category, catcounts in catdict.items():
        # category, catcounts = list(catdict.items())[0]
        # Get total for category
        total = sum(catcounts.values())

        # Work through each read
        for name, count in catcounts.items():
            # Get proportion for this read and add to dict
            if metric == 'p':
                counts[name].add(count / total)
            elif metric == 'n':
                counts[name].add(count)
            else:
                sys.exit("Error: unknown metric passed to score_categories")

    return counts


def reject(counts, threshold, anyfail):
    # counts, threshold = counts[i], t
    # Work through sets of scores
    # name = "uniq13;size=257"
    # count = counts[name]
    out = []
    for name, count in counts.items():
        if ((anyfail and min(count) < threshold)
                or (not anyfail and all(c < threshold for c in count))):
            out.append(name)
    return out


def multicategory(countdict, otherdicts):
    """Find all combinations of sets in category dictionaries.
    First dictionary must have counts of incidences, i.e {x : {a : 1, b : 2}},
    the rest are passed in a tuple of length >= 1, each can have counts or 
    just sets, i.e. {y : {a , b}}
    """
    # countdict, otherdicts = [data[specdict['terms'][0]], parttermsdata]
    # Check otherdicts is a tuple or error

    if type(otherdicts) != tuple:
        otherdicts = (otherdicts,)

    # Initialise master dict with count dict
    multidict = countdict

    # Work through each further dict combining with master
    for newdict in otherdicts:
        # newdict = otherdicts[0]
        # Create new multdict
        newmulti = dict()

        # Work through categories in master
        for mcat, mcounts in multidict.items():
            # mcat, mcounts = list(multidict.items())[0]
            # Work through categories in new dict
            for ncat, nvalue in newdict.items():
                # ncat, nvalue = list(newdict.items())[0]
                # Check if avalue is set of names or dictionary of name:counts
                nvalueset = nvalue
                if type(nvalue) == "dict":
                    nvalueset = nvalue.keys()

                # Find all names shared between current categories
                nnames = set(mcounts.keys()).intersection(nvalueset)

                # Skip combination if no names shared
                if len(nnames) == 0:
                    continue

                # Create the combined category name
                ncat = mcat + ncat

                # Create a combined subdictionary
                newmulti[ncat] = dict()

                # Add the minimum value of the two counts for each shared name
                if type(nvalue) == "dict":
                    newmulti[ncat] = {nn: min(mcounts[nn], nvalue[nn]) for nn in nnames}
                else:
                    newmulti[ncat] = {nn: mcounts[nn] for nn in nnames}

        # Overwrite old multi_dict
        multidict = newmulti

    return multidict


def write_count_dict(countdict, asvs, path):

    # TODO: Check the count dict is of the right structure

    # Order asvs name
    asvs_sort = sorted(asvs)

    # Write to file
    with open(path, 'w') as o:

        # Write header
        o.write(",".join([""] + asvs_sort) + "\n")

        # Write data lines
        for lib, counts in countdict.items():
            values = []
            for asv in asvs_sort:
                values.append(str(counts[asv] if asv in counts else 0))
            o.write(",".join([lib] + values) + "\n")


def resolve_ranges(range_string):
    # Split into segments
    ranges = re.split(",", range_string)
    # Separate parts of any range specifications
    ranges = [re.split("[-/]", r) for r in ranges]
    # TODO: add check for range specifications being complete
    # Convert any range specifications to full list
    rangeout = []
    for r in ranges:
        if len(r) > 1:
            rangeout.append(numpy.linspace(float(r[0]), float(r[1]), int(r[2])))
        else:
            rangeout.append([float(r[0])])
    # Collapse list of lists
    return list(itertools.chain(*rangeout))


def resolve_spec(spec):
    values = re.sub(r'[\[\]]', '', spec).split(';')

    # Set up error
    err = f"Error, malformed specification in {spec}"

    # Check there are 3 values
    if len(values) != 3:
        sys.exit(f"{err}: each specification term should have three semicolon-separated entries, "
                 f"with terms enclosed in square brackets")

    name = f"{values[0]}_{values[1]}"
    # Split terms into parts
    values[0] = re.split(r'[|+]', values[0])

    # Check terms are formed correctly
    if not values[0][0] in ['library', 'total']:
        sys.exit(f"{err}: first term must be \'library\' or \'total\'")
    if len(values[0]) > 1 and not all(t in ['taxon', 'clade'] for t in values[0][1:]):
        sys.exit(f"{err}: secondary term(s) must be \'clade\' or \'taxon\'")

    # Check metric is n or p
    if not values[1] in ['n', 'p']:
        sys.exit(f"{err}: metric should be \'n\' or \'p\'")

    # Expand thresholds to list
    values[2] = resolve_ranges(values[2])

    return [name] + values


def parse_specs(args, null=float('nan')):
    """Parse through the specification to expand terms and thresholds"""
    # null = 0

    if args.mode == 'find':
        # Open specifications file and parse lines into specifications
        fh = open(args.specification, 'r')
        spectext = ''
        for line in fh:
            # line = fh.readline()
            if re.match(r'^\s*#|^$', line):
                continue
            else:
                spectext += line.strip()
        fh.close()
    elif args.mode == 'dump':
        spectext = '*'.join([re.sub(r'\'', '', s) for s in args.specification])
    elif args.mode == 'filter-adaptive':
        # Dummy specs for adaptive mode
        return {'name': []}, set(), 0, 0, []
    else:
        sys.exit(f"Unrecognised args.mode: {args.mode}")

    # Check
    if not args.taxgroups and "taxon" in spectext:
        sys.stderr.write("Warning: 'taxon' specified in binning strategy but no taxon file supplied. "
                         "Ignoring schemes including 'taxon'.\n")

    # Clean up and parse for additive and multiplicative specs
    spectext = re.sub(' ', '', spectext)
    specslist = [s.split('*') for s in spectext.split('+')]

    # Find unique terms and extract thresholds
    specdict = dict()
    termvals = []
    termset = set()
    for specl in specslist:
        out = dict()
        
        # Check for taxon requirement in this additive term
        valid_term = True
        resolved_specs = []
        for spec in specl:
            res = resolve_spec(spec)
            if not args.taxgroups and 'taxon' in res[1]:
                valid_term = False
                break
            resolved_specs.append(res)

        if not valid_term:
            continue

        for res in resolved_specs:
            termset.update(res[1])
            specdict[res[0]] = res[1:3]
            out[res[0]] = res[3]
        termvals.append(out)
    nterm = len(termset)

    specs = {'name': list(specdict.keys())}
    specs.update({k: v for k, v in zip(['spec', 'metric'], map(list, zip(*specdict.values())))})

    # Build all threshold combinations
    threshlists = []
    # termdetails = []
    nthresh = 0
    for tval in termvals:
        # tval = termvals[6]
        specl = []
        threshl = []
        for spec in specs['name']:
            threshl.append(tval[spec] if spec in tval else [null])
            if spec in tval:
                specl.append(spec)
        num = reduce(lambda x, y: x * y, [len(t) for t in threshl], 1)
        nthresh += num
        # termdetails.append([specl, n])
        threshlists.append([specl, threshl])

    if args.mode == 'dump' and nthresh > 1:
        sys.exit("Error: mode 'dump' allows only a single threshold set. To run wih multiple "
                 "threshold sets, use mode 'find'\n")

    # Construct generator
    # termgen = (l for t, c in termdetails for l in [t] * c)
    # thresholds = (t for tl in threshlists for t in itertools.product(*tl))
    def threshold_gen(thrlis):
        n = 0
        for tn, tl in thrlis:
            for t in itertools.product(*tl):
                yield tuple([n, tn]) + t
                n += 1

    thresholds = threshold_gen(threshlists)

    return specs, termset, nterm, nthresh, thresholds


def get_validated(raw, args, basepath, totalcounts):
    # filname = infilename
    # Set up path
    controlpath = basepath + "_control.txt"

    # Check if can resume
    if os.path.exists(controlpath) and not args.overwrite:
        sys.stdout.write(f"Resuming with validated ASV lists found in {controlpath}.\n")
        target, nontarget = [], []
        with open(controlpath, 'r') as i:
            for line in i:
                cat, asv = line.rstrip().split('\t')
                if cat == 'refpass':
                    target.append(asv)
                else:
                    nontarget.append(asv)
        sys.stdout.write(f"Read {len(target)} validated authentic ASVs and {len(nontarget)} "
                         f"validated non-authentic ASVs from {controlpath}.\n")
        return set(target), set(nontarget)

    # Create length-based control list
    sys.stdout.write("Identifying validated non-authentic ASVs based on length...")

    nontargetlength = filterlength.check_length_multi(raw['asvs'],
                                                      [args.minimumlength,
                                                       args.maximumlength],
                                                      args, fail=True)

    sys.stdout.write(f"found {len(nontargetlength)} candidates\n")

    # Create translation based control list
    sys.stdout.write("Identifying validated non-authentic ASVs based on translation...")

    try:
        nontargettrans = filtertranslate.check_stops_multi(raw['asvs'], args, fail=True)
    except BaseException as e:
        sys.stdout.write(f" Warning: Translation filter failed ({e}), skipping...\n")
        nontargettrans = []

    sys.stdout.write(f"found {len(nontargettrans)} candidates\n")

    # Finalise nontargets
    nontarget = set(nontargetlength + nontargettrans)
    overlap = set(nontargetlength).intersection(set(nontargettrans))
    if len(nontarget) > 0:
        sys.stdout.write(f"Validated a total of {len(nontarget)} unique ASVs as non-authentic "
                         f"({len(overlap)} were found in both sets)\n")
    else:
        sys.exit("Error: no non-authentic ASVs could be found. Prior data filtering may have been "
                 "too stringent.")

    # Create reference based control list

    wd = filterreference.make_temp_bbmapwd(args.output, "bbmap")
    allcandidates = []

    sys.stdout.write(f"Identifying validated authentic ASVs based on {args.references}...")
    sys.stdout.flush()

    # mp, ml = [args.refmatchpercent, args.refmatchlength]
   
    candidates = filterreference.refmatch_BBMap(raw['path'], wd, args.refmatchlength, args.threads, args.references, totalcounts, args)

    allcandidates.extend(candidates)
    
    if not args.keeptemporaryfiles:
        shutil.rmtree(wd)

    refmatch = set(allcandidates)
    target = refmatch - nontarget

    # Finalise targets
    if len(target) > 0:
        sys.stdout.write(f"Validated a total of {len(target)} unique ASVs as authentic "
                         f"({len(refmatch)} unique matched to references, "
                         f"{len(refmatch - target)} rejected due to inclusion in non-authentic "
                         f"set)\n")
    else:
        err = "Error: no authentic ASVs found"
        if len(refmatch) > 0:
            err = (f"{err}, although {len(refmatch)} ASVs matched to reference set. Length "
                   f"thresholds may be too stringent and find too many non-authentic ASVs, or "
                   f"reference dataset is not sufficiently curated")
        else:
            err = (f"{err}. Check the reference file and/or database is correct and consider "
                   f"adjusting the matching thresholds (carefully!)")
        sys.exit(err)

    # Output file with details of targets/non-targets
    with open(controlpath, 'w') as o:
        for cat, asv in zip(["lengthfail", "stopfail", "refpass"],
                            [nontargetlength, nontargettrans, refmatch]):
            for a in asv:
                o.write(f"{cat}\t{a}\n")

    return target, nontarget


def counts_from_spec(specs, data):
    """Count categories for each individual specification in a specification set and return a list
    of the category counts for each specification"""

    counts = dict()

    # Work through specifications
    names = specs.get('name', [])
    sps = specs.get('spec', [])
    metrics = specs.get('metric', [])

    for name, spec, metric in zip(names, sps, metrics):
        # name, spec, metric = list(zip(*specs.values()))[1]
        # Check if partitioned and do counting
        if len(spec) == 1:
            counts[name] = count_categories(data[spec[0]], metric)
        else:
            partspecdata = tuple(data[s] for s in spec[1:])
            multicat = multicategory(data[spec[0]], partspecdata)
            counts[name] = count_categories(multicat, metric)

    return counts


def estimate_true_values(asvs, retained_asvs, retained_target, target, retained_nontarget,
                         nontarget):
    try:
        true_target = ((retained_asvs - asvs * (retained_nontarget / nontarget)
                        ) / (retained_target / target - retained_nontarget / nontarget))
        true_nontarget = asvs - true_target
        true_retained_target = true_target * (retained_target / target)
        true_retained_nontarget = ((retained_nontarget / nontarget) * (asvs - true_target))
        proportion_retained_target = ( true_retained_target / retained_asvs )
        proportion_retained_nontarget = ( true_retained_nontarget / retained_asvs ) 
        out = [true_target, true_nontarget, true_retained_target, true_retained_nontarget, proportion_retained_target, proportion_retained_nontarget]
    except:
        out = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA']
    return out


def calc_scores(truepositive, falsepositive, truenegative, falsenegative):
    def _divide(num, den):
        return num / den if den else 0

    accuracy = _divide((truepositive + truenegative),
                       (truepositive + truenegative + falsepositive + falsenegative))
    precision = _divide(truepositive, (truepositive + falsepositive))
    recall = _divide(truepositive, (truepositive + falsenegative))

    return [accuracy, precision, recall]


def apply_reject(threshnames, thresholds, counts, anyfail):
    rejects = []
    for term, thresh in zip(threshnames, thresholds):
        # term, thresh = list(zip(threshnames, thresholds))[3]
        rejects.extend(reject(counts[term], thresh, anyfail))
    return set(rejects)


def calc_stats(rejects, asvs, target, nontarget):
    """For a given set of category counts and a given set of thresholds, counts retention,
    calculates scores and estimates statistics. Counts and thresholds should be in two lists of
    equal lengths, where the nth item of the thresholds should be the threshold count for the
    nth counts"""

    # Find all rejected_asvs across counts
    actualrejects = rejects.union(nontarget) - target

    # Calculate summary statistics
    stats = [len(asvs),  # 0: asvs_total
             len(target),  # 1: targets_total
             len(nontarget),  # 2: nontargets_total
             len(asvs) - len(rejects)]  # 3: asvsprelim_retained_n
    stats += [stats[3] / stats[0],  # 4: asvsprelim_retained_p
              len(rejects),  # 5: asvsprelim_rejected_n
              len(rejects) / stats[0],  # 6: asvsprelim_rejected_p
              len(target - rejects)]  # 7: targets_retained_n
    stats += [stats[7] / stats[1],  # 8: targets_retained_p
              stats[1] - stats[7]]  # 9: targets_rejected_n
    stats += [stats[9] / stats[1],  # 10: targets_rejected_p
              len(nontarget - rejects)]  # 11: nontargets_retained_n
    stats += [stats[11] / stats[2],  # 12: nontargets_retained_p
              stats[2] - stats[11]]  # 13: nontargets_rejected_n
    stats += [stats[13] / stats[2],  # 14: nontargets_rejected_p
              len(asvs) - len(actualrejects)]  # 15: asvsactual_retained_n
    stats += [stats[15] / stats[0],  # 16: asvsactual_retained_p
              len(actualrejects),  # 17: asvsactual_rejected_n
              len(actualrejects) / stats[0]]  # 18: asvsactual_rejected_p

    # True positives = targets_retained_n = 7
    # False positives = nontargets_retained_n = 11
    # True negatives = nontargets_rejected_n = 13
    # False negatives = targets_rejected_n = 9

    # Calculate score from truepositives, falsepositives, truenegatives and
    # falsenegatives
    scores = calc_scores(stats[7], stats[11], stats[13], stats[9])

    # Calculate estimates of 0: total input true targets, 1: total input true
    # nontargets, 2: total output true targets, 3: total output true nontargets
    estimates = estimate_true_values(stats[0], stats[3], stats[7], stats[1],
                                     stats[11], stats[2])

    # Store list of rejects or retains, whichever is shorter
    if len(actualrejects) < 0.5 * len(asvs):
        store = ('reject', actualrejects)
    else:
        store = ('retain', asvs - actualrejects)

    # score, stats, estimates, hash, store
    rejects = tuple(sorted(actualrejects))
    return scores, stats + estimates, hash(rejects), store


def assess_numts(threshnames, counts, anyfail, asvs, target, nontarget, queues, threshvals):
    # threshnames, anyfail, asvs = specs['name'], args.anyfail, set(raw['asvs'].keys())
    # threshvals = next(thresholds)

    # Parse link data
    i, terms, threshvals = threshvals[0], threshvals[1], threshvals[2:]

    # Run rejection and statistics
    rejects = apply_reject(threshnames, threshvals, counts, anyfail)
    scores, stats, ohash, store = calc_stats(rejects, asvs, target, nontarget)

    # Write to terminal and files
    if queues:
        prinq, termq = queues
        prinq.put((i, terms, list(threshvals) + scores + stats, ohash, store))
        termq.put(ohash)

    return [i] + scores


def write_stats_and_cache(specs, outdir, prinq):
    sh = open(os.path.join(outdir, "results.csv"), 'w')
    ch = gzip.open(os.path.join(outdir, "resultcache"), 'wt')
    hh = open(os.path.join(outdir, "hashcache"), 'w')

    # Write header
    if specs.get('otu_mode'):
        ch.write("# otu_mode=True\n")
    
    head = ("accuracy_score precision_score recall_score "
            "asvs_total verifiedauthentic_total_observed "
            "verifiednonauthentic_total_observed "
            "asvsprelim_retained_n asvsprelim_retained_p "
            "asvsprelim_rejected_n asvsprelim_rejected_p "
            "verifiedauthentic_retained_n verifiedauthentic_retained_p "
            "verifiedauthentic_rejected_n verifiedauthentic_rejected_p "
            "verifiednonauthentic_retained_n verifiednonauthentic_retained_p "
            "verifiednonauthentic_rejected_n verifiednonauthentic_rejected_p "
            "asvsactual_retained_n asvsactual_retained_p "
            "asvsactual_rejected_n asvsactual_rejected_p "
            "authentic_total_estimate nonauthentic_total_estimate "
            "authentic_retained_estimate_n nonauthentic_retained_estimate_n "
            "authentic_retained_estimate_p nonauthentic_retained_estimate_p "
            "rejects_hash").split(" ")

    sh.write(",".join(["resultindex", "term"]
                      + [s + "_threshold" for s in specs['name']]
                      + head) + '\n')

    while 1:
        feed = prinq.get()
        if feed is None:
            break
        i, term, statl, ohash, store = feed
        statl = [i, '*'.join(term)] + statl + [ohash]
        sh.write(",".join(str(v) for v in statl) + '\n')
        sh.flush()
        ch.write('\t'.join([str(i), store[0]] + list(store[1])) + '\n')
        ch.flush()
        hh.write(f"{i}\t{ohash}\n")
        hh.flush()

    sh.close()
    ch.close()
    hh.close()


def write_terminal(tot, termq):
    # Set up to print
    start = time.perf_counter()
    done = 0
    remain = "unknown time"
    uniqouts = set()
    # Print
    while 1:
        sys.stdout.write("\r%s\r" % (' ' * 70))
        sys.stdout.write(f"\rAssessed {done} of {tot} total threshold combinations, {remain} "
                         f"remaining")
        sys.stdout.flush()
        queueitem = termq.get()
        if queueitem is None:
            break
        uniqouts.add(queueitem)
        done += 1
        now = time.perf_counter()
        elapsed = now - start
        remain = round((elapsed / done) * (tot - done))
        remain = "approx " + str(datetime.timedelta(seconds=remain))

    now = time.perf_counter()
    elapsed = now - start
    elapsedper = datetime.timedelta(seconds=elapsed / done)
    elapsed = datetime.timedelta(seconds=round(elapsed))
    sys.stdout.write(f"\nFinished assessing {tot} threshold combinations in {elapsed}, "
                     f"{elapsedper} per threshold combination.\n"
                     f"Generated a total of {len(uniqouts)} unique ASV output compositions\n")
    sys.stdout.flush()


def start_writers(pool, manager, specs, outdir, nthresh):
    prinq = manager.Queue()
    printwatch = pool.apply_async(partial(write_stats_and_cache, specs, outdir), (prinq,))

    termq = manager.Queue()
    termwatch = pool.apply_async(partial(write_terminal, nthresh), (termq,))

    return (prinq, termq), (printwatch, termwatch)


def get_reject_from_store(asvs, n, store):
    # store = rsstore
    if store[0] == 'retain':
        return asvs - store[1]
    elif store[0] == 'reject':
        return store[1]
    else:
        sys.exit(f"store type \'{store[0]}\' for resultset {n} is not \'retain\' or \'reject\'")


def find_best_score(scores, scoretype):
    # scoretype = args.scoremetric
    scoreloc = {'accuracy': 0,
                'precision': 1,
                'recall': 2}

    scorelist = [[s[0], s[1 + scoreloc[scoretype]]] for s in scores]
    scorelist = sorted(scorelist, key=lambda x: int(x[1]))[::-1]
    maxscore = scorelist[0][1]
    nmax = len([s for s in scorelist if s[1] == maxscore])
    sys.stdout.write(f"Maxmimum {scoretype} score of {round(maxscore, 4)} achieved by {nmax} "
                     f"threshold sets\n")

    return scorelist


def parse_hashcache(path):
    hashdict = defaultdict(set)
    fh = open(path, 'r')
    for line in fh:
        rs, h = line.strip().split('\t')
        hashdict[h].add(int(rs))
    fh.close()
    return hashdict


def generate_resultsets(genval, scoresort, path):
    # genval, path = args.generateASVresults, hashcachepath
    uniqscores = sorted(list(set([s[1] for s in scoresort])))[::-1]
    i = 0
    p = genval
    if type(genval) is float:
        i = round(genval * len(uniqscores)) - 1
        p = f"{round(genval * 10)}% of"
    scores = uniqscores[:i + 1]
    minscore, maxscore = min(scores), max(scores)

    resultsets = [s[0] for s in scoresort if s[1] >= minscore]

    hashcache = parse_hashcache(path)
    users = []
    for h, rs in hashcache.items():
        # h, rs = list(hashcache.items())[0]
        inrs = [r for r in rs if r in resultsets]
        if len(inrs) > 0:
            users.append(inrs[0])

    sys.stdout.write(f"Found {len(uniqscores)} unique scores, outputting the results for the top "
                     f"{p} score(s) (score(s) {minscore} to {maxscore}), corresponding to "
                     f"{len(resultsets)} threshold combinations and {len(users)} unique ASV "
                     "output compositions\n")

    return users


def write_retained_asvs(infile, outfile, rejects):
    with open(infile, 'r') as infa, open(f"{outfile}", 'w') as outfa:
        for head, seq in SimpleFastaParser(infa):
            if head not in rejects:
                outfa.write(f">{head}\n{seq}\n")


def enforce_validation_fasta(fasta_path, source_path, target, nontarget):
    """Post-process a FASTA: add missing authentic ASVs, remove non-authentic ASVs."""
    # Read current output
    current = {}
    with open(fasta_path, 'r') as f:
        for head, seq in SimpleFastaParser(f):
            current[head] = seq

    # Remove non-authentics
    for name in nontarget:
        current.pop(name, None)

    # Add missing authentics from source
    with open(source_path, 'r') as f:
        for head, seq in SimpleFastaParser(f):
            if head in target and head not in current:
                current[head] = seq

    # Rewrite
    with open(fasta_path, 'w') as f:
        for head, seq in current.items():
            f.write(f">{head}\n{seq}\n")


def load_validation_from_control(controlpath):
    """Load target and nontarget sets from a control file."""
    target, nontarget = [], []
    with open(controlpath, 'r') as f:
        for line in f:
            cat, asv = line.rstrip().split('\t')
            if cat == 'refpass':
                target.append(asv)
            else:
                nontarget.append(asv)
    return set(target), set(nontarget)


def load_validation_from_otu_summary(summarypath):
    """Load target and nontarget OTU sets from an otu_summary.csv file."""
    target, nontarget = set(), set()
    with open(summarypath, 'r') as f:
        header = f.readline()
        for line in f:
            parts = line.rstrip().split(',')
            otu, otu_status = parts[1], parts[3]
            if otu_status == 'Authentic':
                target.add(otu)
            elif otu_status == 'Non-Authentic':
                nontarget.add(otu)
    return target, nontarget


def make_resultset_paths(name, resultsets):
    if len(resultsets) > 1:
        paths = {rs: f"{name}_resultset{rs}.fasta" for rs in resultsets}
    else:
        if os.path.splitext(name)[1]:
            paths = {resultsets[0]: name}
        else:
            paths = {resultsets[0]: f"{name}.fasta"}
    return paths


def write_resultset_asvs(asvs, resultsets, cachepath, infile, name, mode):

    # Read in the relevant store lines
    store = parse_resultcache(cachepath, asvs, resultsets)
    if mode == 'dump':
        sys.stdout.write(f"Read {len(store)} cached results.\n")

    sys.stdout.write(f"Writing fastas for {len(resultsets)} results...")

    paths = make_resultset_paths(name, resultsets)

    for rs, rsstore in store:
        # rs, rsstore = store[0]
        rejects = get_reject_from_store(asvs, rs, rsstore)
        write_retained_asvs(infile, paths[rs], rejects)
    sys.stdout.write("done.\n")


def parse_resultcache(path, asvs, resultsets=None):
    # path, asvs = [args.resultcache, set(raw['asvs'].keys())]
    fh = gzip.open(path, 'rt')
    store = []
    for i, line in enumerate(fh):
        # i, line = next(enumerate(fh))
        if line.startswith('#'):
            continue
        vals = line.strip().split('\t')
        err = f"Error: {path} line {i + 1}"
        if len(vals) < 2:
            sys.exit(f"{err} does not have enough items")
        setn = vals.pop(0)
        try:
            setn = int(setn)
        except:
            sys.exit(f"{err} starting value '{setn}' is not an integer")
        if resultsets and setn not in resultsets:
            continue
        action = vals.pop(0)
        if action not in {'reject', 'retain'}:
            sys.exit(f"{err} second value '{action}' is not 'reject' or 'retain'")
        missing = [a for a in vals if a not in asvs]
        if len(missing) > 0:
            sys.exit(f"{err} ASVs {', '.join(missing)} are not present or identifiable in the "
                     f"supplied ASV file")
        store.append([setn, (action, set(vals))])
    fh.close()

    return store


def check_cache_otu_mode(path):
    """Check if the cache was created in OTU mode."""
    mode = False
    try:
        with gzip.open(path, 'rt') as fh:
            line = fh.readline()
            if line.startswith('# otu_mode=True'):
                mode = True
    except Exception:
        pass
    return mode


def adaptive_filter(asvs, librarycounts, target, nontarget, percentile, criteria, output_csv, output_fasta, output_summary):
    """
    Perform per-sample filtering based on the distribution of non-authentic ASVs.
    
    Args:
        asvs: dictionary of ASV identifiers (metadata/sequences)
        librarycounts: dictionary mapping library names to dictionaries of {asv: count}
        target: set of authentic ASVs
        nontarget: set of non-authentic ASVs
        percentile: percentile to target (e.g. 0.95)
        criteria: 'verified_removed' or 'estimated_removed'
        output_csv: path to write filtered abundance table
        output_fasta: path to write valid ASV sequences
        output_summary: path to write summary statistics table
    """
    
    valid_thresholds = []
    sample_thresholds = {}
    
    # 1. Calculate per-sample thresholds
    skipped_samples = [] # Samples where we couldn't calculate a threshold
    
    for lib, counts in librarycounts.items():
        # Get counts for this library
        
        # Identify non-targets in this library
        lib_nontargets = [count for asv, count in counts.items() if asv in nontarget]
        
        threshold = None
        
        if criteria == 'verified_removed':
            if len(lib_nontargets) > 0:
                # Calculate percentile of verified non-authentic abundances
                # If we want to remove X% of non-authentic ASVs, we find the Xth percentile value.
                # All ASVs with count < threshold will be removed (check numpy.percentile behavior)
                # numpy.percentile with 95 means 95% of data is below result.
                # So if we want to remove 95%, we use percentile * 100.
                threshold = numpy.percentile(lib_nontargets, percentile * 100)
            else:
                 # No verified non-authentic ASVs to base threshold on
                 pass
                 
        elif criteria == 'estimated_removed':
             # Get unique counts to speed up the process
             unique_counts = sorted(list(set(counts.values())))
             
             best_diff = float('inf')
             best_t = None
             
             # Prepare total stats for this library to save re-calculation 
             lib_total_asvs = len(counts)
             
             if lib_total_asvs < 10: # Minimum data check
                  pass
             else:
                 # Check if we have enough authentic/non-authentic markers
                 lib_target_set = set(counts.keys()).intersection(target)
                 lib_nontarget_set = set(counts.keys()).intersection(nontarget)
                 
                 if len(lib_target_set) > 0 and len(lib_nontarget_set) > 0:
                     for t in unique_counts:
                         # Simulate filtering at threshold t
                         # Retained: >= t
                         # Rejected: < t
                         
                         retained_keys = [k for k, v in counts.items() if v >= t]
                         # rejected_keys = [k for k, v in counts.items() if v < t]
                         
                         retained_target_n = len(set(retained_keys).intersection(target))
                         retained_nontarget_n = len(set(retained_keys).intersection(nontarget))
                         
                         retained_asvs_n = len(retained_keys)
                         
                         # Global stats for this library needed for estimation
                         target_n = len(lib_target_set)
                         nontarget_n = len(lib_nontarget_set)
                         
                         # Run estimation
                         est = estimate_true_values(lib_total_asvs, retained_asvs_n, retained_target_n, target_n, retained_nontarget_n, nontarget_n)
                         
                         true_nontarget = est[1]
                         true_retained_nontarget = est[3]
                         
                         if true_nontarget != 'NA' and true_nontarget > 0:
                             # Calculate removed proportion
                             removed_prop = 1.0 - (true_retained_nontarget / true_nontarget)
                             
                             diff = abs(removed_prop - percentile)
                             if diff < best_diff:
                                 best_diff = diff
                                 best_t = t
                     
                     threshold = best_t

        if threshold is not None:
             valid_thresholds.append(threshold)
             sample_thresholds[lib] = threshold
        else:
             skipped_samples.append(lib)
             
             
    # 2. Handle Fallback
    fallback_threshold = 0
    if valid_thresholds:
        fallback_threshold = numpy.mean(valid_thresholds)
        sys.stdout.write(f"Calculated valid thresholds for {len(valid_thresholds)} samples. Mean: {fallback_threshold:.2f}\n")
    else:
        sys.stdout.write("Warning: No valid thresholds calculated using specified criteria. Defaulting to 1 (keep all with count >= 1).\n")
        fallback_threshold = 1
        
    for lib in skipped_samples:
        sample_thresholds[lib] = fallback_threshold
        
    # 3. Apply Filtering and Collect Summary Stats
    filtered_library_counts = {}
    retained_asvs_set = set()
    summary_data = [] # List of tuples/dicts to write to summary CSV
    
    for lib, counts in librarycounts.items():
        thresh = sample_thresholds[lib]
        
        # Pre-filter stats
        asvs_present = set(counts.keys())
        total_asvs_pre = len(asvs_present)
        va_present = asvs_present.intersection(target)
        vna_present = asvs_present.intersection(nontarget)
        va_count_pre = len(va_present)
        vna_count_pre = len(vna_present)
        
        # Apply filter (threshold only; enforcement is applied as post-processing)
        new_counts = {asv: c for asv, c in counts.items() if c >= thresh}
        filtered_library_counts[lib] = new_counts
        
        # Post-filter stats
        retained_keys = set(new_counts.keys())
        total_asvs_post = len(retained_keys)
        # va_retained = retained_keys.intersection(va_present) # Equivalent to intersection(target)
        # vna_retained = retained_keys.intersection(vna_present) # Equivalent to intersection(nontarget)
        va_count_post = len(retained_keys.intersection(target))
        vna_count_post = len(retained_keys.intersection(nontarget))
        
        retained_asvs_set.update(retained_keys)
        
        summary_data.append({
            "Sample": lib,
            "Threshold": round(thresh, 2),
            "Total_ASVs_Pre": total_asvs_pre,
            "Verified_Authentic_Pre": va_count_pre,
            "Verified_NonAuthentic_Pre": vna_count_pre,
            "Total_ASVs_Post": total_asvs_post,
            "Verified_Authentic_Post": va_count_post,
            "Verified_NonAuthentic_Post": vna_count_post
        })
        
    # 4. Write CSV
    all_retained_asvs = sorted(list(retained_asvs_set))
    
    with open(output_csv, 'w') as f:
         # Header
         f.write("Sample,Threshold," + ",".join(all_retained_asvs) + "\n")
         
         for lib, counts in filtered_library_counts.items():
             thresh = sample_thresholds.get(lib, fallback_threshold)
             row = [lib, str(round(thresh, 2))]
             for asv in all_retained_asvs:
                 row.append(str(counts.get(asv, 0)))
             f.write(",".join(row) + "\n")
             
    # 5. Write FASTA
    # In metamate.py: raw['asvs'] = SeqIO.to_dict(SeqIO.parse(args.asvs, "fasta"))
    
    with open(output_fasta, 'w') as f:
        for asv_name in all_retained_asvs:
            if asv_name in asvs:
                 record = asvs[asv_name]
                 # If it is Bio.SeqRecord
                 f.write(f">{asv_name}\n{str(record.seq)}\n") # simplified
            else:
                 pass
                 
    # 6. Write Summary Table
    with open(output_summary, 'w') as f:
        headers = ["Sample", "Threshold", "Total_ASVs_Pre", "Verified_Authentic_Pre", "Verified_NonAuthentic_Pre", 
                   "Total_ASVs_Post", "Verified_Authentic_Post", "Verified_NonAuthentic_Post"]
        f.write(",".join(headers) + "\n")
        
        # Sort summary data by sample name for niceness
        summary_data.sort(key=lambda x: x["Sample"])
        
        for row in summary_data:
            line_vals = [str(row[h]) for h in headers]
            f.write(",".join(line_vals) + "\n")