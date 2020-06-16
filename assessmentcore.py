#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core functions for designating numts and optimising parameters"""

# Imports

import categorycounting
import filterlength
import filtertranslate
import filterreference

from Bio.SeqIO.FastaIO import SimpleFastaParser
import re
import numpy
import itertools
import os
import sys

# Global variables

# Inputs:
    # Bad set, as one or more sets of names
    #bad = {'A', 'B', 'C'}

    # Good set, as one or more sets of names
    #good = {'X', 'Y', 'Z'}

    # Counts dictionaries (for libraries or totals) dict[library][name] : count or dict["total"][name] : count
    #library_counts = { "l1" : { "A" : 3, "B" : 2, "D" : 5}}
    #total_counts = {"total" : {"A" : 3, "B" : 2, "D" : 5}}

    # Category dictionaries (for clades or taxa) dict[category] : {names}
    #clades = {"c1" : {"A", "B"}, "c2" : {"D"}}

    #Specifications, as follows:

    # spec = a single specification for a count_categories run, i.e. one each of terms, metric and thresholds
    # spec set = a set of specifications for a single run (i.e. the dictionary below)

    # specs[0][terms] : "total"
    # specs[0][metric] : "n"
    # specs[0][thresholds] : "1-5,1"
    # specs[1][terms] : "library"
    # specs[1][metric] : "p"
    # specs[1][thresholds] : "0.1-0.5,0.05"

    # terms should be specified as "a", "a | b", " a | b + c" where a is the count_dict to use and b (and c) are the category dicts

# Class definitions

# Function definitions

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
            rangeout.append(numpy.linspace(float(r[0]),
                                           float(r[1]),
                                           int(r[2])))
        else:
            rangeout.append([float(r[0])])
    # Collapse list of lists
    return(list(itertools.chain(*rangeout)))

def parse_spec(spec):

    values = re.sub("[\[\]]", '', spec).split(';')
    name = f"{values[0]}_{values[1]}"
    # Set up error
    err = f"Error, malformed specification in {spec}"

    # Check there are 3 values
    if len(values) != 3:
        sys.exit(f"{err}: specification line should have three "
                 "tab-separated entries")

    # Split terms into parts
    values[0]  = re.split("[\|\+]", values[0])

    # Check terms are formed correctly
    if not values[0][0] in ['library', 'total']:
        sys.exit(f"{err}: first term must be \'library\' or \'total\'")
    if (len(values[0]) > 1
        and not all(t in ['taxon', 'clade'] for t in values[0][1:])):
        sys.exit(f"{err}: secondary term(s) must be \'clade\' or "
                 "\'taxon\'")

    # Check metric is n or p
    if(not values[1] in ['n', 'p']):
        sys.exit(f"{err}: metric should be \'n\' or \'p\'")

    # Expand thresholds to list
    values[2] = resolve_ranges(values[2])

    return([name] + values)

def parse_specfile(args, null = float('nan')):
    "Parse through the specification to expand terms and thresholds"

    if args.action == 'find':
        # Open specifications file and parse lines into specifications
        fh = open(args.specification, 'r')
        spectext = ''
        for l in fh:
            #l = fh.readline()
            if re.match("^\s*#|^$", l):
                continue
            else:
                spectext += l.strip()
        fh.close()
    elif args.action == 'dump':
        spectext = '*'.join([re.sub("\'", '', s) for s in args.specification])
    
    # Clean up and parse for additive and multiplicative specs
    spectext = re.sub(' ', '', spectext)
    specslist = [s.split('*') for s in spectext.split('+')]

    # Find unique terms and extract thresholds
    specdict = dict()
    specthresholds = []
    for specl in specslist:
        out = dict()
        for spec in specl:
            l = parse_spec(spec)
            specdict[l[0]] = {'terms': l[1],
                              'metric': l[2]}
            out[l[0]] = l[3]
        specthresholds.append(out)
    # Build all threshold combinations
    threshcombos = []
    termorder = list(specdict.keys())
    for specs in specthresholds:
        #specs = specthresholds[2]
        thresholds = []
        for term in termorder:
            thresholds.append(specs[term] if term in specs else [null])
        threshcombos.extend(itertools.product(*thresholds))
    # Check
    if( not args.taxgroups and any("taxon" in e for e in spectext)):
        sys.exit("Error, taxon specified as a binning strategy but no taxon "
                 "file supplied")
    #Output specifications in table
#    if args.action == 'find':
#        filename = os.path.splitext(os.path.basename(args.asvs))[0]
#        path = os.path.join(args.outputdirectory,
#                            f"{filename}_specifications.csv")
#        with open(path, 'w') as o:
#            o.write(",".join(specdict.keys()) + "\n")
#            for thresh in threshcombos:
#                o.write(",".join([str(t) for t in thresh]) + "\n")
#
    return(specdict, termorder, threshcombos)

def get_validated(raw, args, filename):
    #filname = infilename
    
    # Create length-based control list
    sys.stdout.write("Identifying control non-target ASVs based on length.\n")

    nontargetlength = filterlength.check_length_multi(raw['asvs'],
                                                      [args.minimumlength,
                                                      args.maximumlength],
                                                      args, fail=True)

    # Create translation based control list
    sys.stdout.write("Identifying control non-target ASVs based on "
                     "translation.\n")

    nontargettrans = filtertranslate.check_stops_multi(raw['asvs'], args,
                                                       fail = True)

    # Finalise nontargets
    nontarget = set(nontargetlength + nontargettrans)

    if len(nontarget) > 0:
        sys.stdout.write(f"Found {len(nontarget)} non-target ASVs: "
                         f"{len(nontargetlength)} based on length variation "
                         f"and {len(nontargettrans)} based on translation.\n")
    else:
        sys.exit("Error: no non-target ASVs could be found. Prior data "
                 "filtering may have been too stringent.")

    # Create reference based control list
    
    
    wd = filterreference.make_temp_blastwd(args.outputdirectory, "blastdb")
    refmatch = set()
    loopvars = zip(['references', 'blast database'], 
                   [args.references, args.blastdb])
    for src, dat in loopvars:
        #src, dat = list(loopvars)[0]
        if dat is None:
            continue
        
        sys.stdout.write(f"Identifying control target ASVs based on {src}...")
        sys.stdout.flush()
        db, mp, ml = [None, None, None]
        
        if src == 'references':
            db = filterreference.make_blastdb(dat, wd) 
            mp, ml = [args.refmatchpercent, args.refmatchlength]
        else:
            db = dat
            mp, ml = [args.dbmatchpercent, args.dbmatchlength]
        
        candidates = filterreference.refmatch_blast(raw['path'], db, wd,
                                                  mp, ml, args.threads)
        sys.stdout.write(f"found {len(candidates)} candidates\n")
        refmatch.update(candidates)
    
    target = refmatch - nontarget
    
    # Finalise targets
    if len(target) > 0:
        sys.stdout.write(f"Found {len(target)} target ASVs: {len(refmatch)} "
                          "out of all ASVs matched to references and/or "
                         f"blast database, {len(refmatch - target)} of these "
                          "rejected due to inclusion in non-target set\n")
    else:
        err = "Error: no target ASVs found"
        if len(refmatch) > 0:
            err = (f"{err}, although {len(refmatch)} ASVs matched to "
                    "reference set. Length thresholds may be too stringent "
                    "and find too many non-target ASVs.")
        else:
            err = (f"{err}. Check the reference file and/or database is "
                   " correct and consider adjusting the matching "
                   " thresholds (carefully!)")
        sys.exit(err)

    # Output file with details of targets/non-targets
    path = os.path.join(args.outputdirectory, f"{filename}_control.txt")
    with open(path, 'w') as o:
        for cat, asv in zip(["lengthfail", "stopfail", "refpass"],
                            [nontargetlength, nontargettrans, refmatch]):
            for a in asv: o.write(f"{cat}\t{a}\n")

    return(target, nontarget)


def counts_from_spec(spec, data):
    """Count categories for each individual specification in a specification
    set and return a list of the category counts for each specification"""
    # spec = specs


    counts = dict()

    # Work through specifications
    for specname, specdict in spec.items():
        #specname, specdict = list(spec.items())[2]
        # Check if partitioned and do counting
        if len(specdict['terms']) == 1 :
            counts[specname] = categorycounting.count_categories(
                                                    data[specdict['terms'][0]],
                                                    specdict['metric'])
        else:
            parttermsdata = tuple( data[t] for t in specdict['terms'][1:] )
            multicat = categorycounting.multicategory(
                                                    data[specdict['terms'][0]],
                                                    parttermsdata)
            counts[specname] = categorycounting.count_categories(multicat,
                                                            specdict['metric'])

    return(counts)

def calc_score(retained_target, target, retained_nontarget, nontarget,
               score_type, weight = 0.5):
    if(score_type == "standardised"):
        return(2 * ( weight * (1 - retained_target/target)
               + (1-weight) * retained_nontarget/nontarget))
    elif(score_type == "unstandardised"):
        return(( weight * (target - retained_target)
                + (1 - weight) * (nontarget - retained_target))
               / ( weight * target + (1 - weight) * nontarget))
    else:
        sys.ext("Error: unknown score_type value passed to calc_score")

def estimate_true_values(asvs, retained_asvs, retained_target, target,
                         retained_nontarget, nontarget):
    try:
        true_target = ((retained_asvs
                        - asvs * (retained_nontarget / nontarget)
                       ) / (retained_target / target
                            - retained_nontarget / nontarget))
        true_nontarget = asvs - true_target
        true_retained_target = true_target * (retained_target / target)
        true_retained_nontarget = ((retained_nontarget / nontarget)
                                   * (asvs - true_target))
        out = [true_target, true_nontarget,
               true_retained_target, true_retained_nontarget]
    except:
        out = ['NA', 'NA', 'NA', 'NA']
    return(out)

def assess_numts(termorder, counts, anyfail,
                 asvs, target, nontarget, scoretype, weight,
                 thresholds):
    rejects = apply_reject(termorder, thresholds, counts, anyfail)
    stats = calc_stats(rejects, asvs, target, nontarget, scoretype, weight)
    return(stats)

def apply_reject(termorder, thresholds, counts, anyfail):
    rejects = []
    for term, thresh in zip(termorder, thresholds):
        #term, thresh = list(zip(termorder, thresholds))[3]
        rejects.extend(categorycounting.reject(counts[term], thresh, anyfail))
    return(set(rejects))

def calc_stats(rejects, asvs, target, nontarget, scoretype, weight):
    #asvs, anyfail, scoretype, thresholds, weight =  [set(raw['asvs'].keys()), args.anyfail, "standardised" , thresholdcombos[0], 0.5]
    """For a given set of category counts and a given set of thresholds, counts
    retention, calculates scores and estimates statistics. Counts and
    thresholds should be in two lists of equal lengths, where the nth item of
    the thresholds should be the threshold count for the nth counts"""

    # TODO: set up multithreading to write to a specs file in parallel to
    # save memory.

    # Find all rejected_asvs across counts
    actualrejects = rejects - target

    stats = [len(asvs),                              #0: asvs_total
             len(target),                            #1: targets_total
             len(nontarget),                         #2: nontargets_total
             len(asvs) - len(rejects)]               #3: asvsprelim_retained_n
    stats += [stats[3] / stats[0],                   #4: asvsprelim_retained_p
              len(rejects),                          #5: asvsprelim_rejected_n
              len(rejects) / stats[0],               #6: asvsprelim_rejected_p
              len(target - rejects)]                 #7: targets_retained_n
    stats += [stats[7] / stats[1],                   #8: targets_retained_p
              stats[1] - stats[7]]                   #9: targets_rejected_n
    stats += [stats[9] / stats[1],                   #10: targets_rejected_p
              len(nontarget - rejects)]              #11: nontargets_retained_n
    stats += [stats[11] / stats[2],                  #12: nontargets_retained_p
              stats[2] - stats[11]]                  #13: nontargets_rejected_n
    stats += [stats[13] / stats[2],                  #14: nontargets_rejected_p
              len(asvs) - len(actualrejects)]        #15: asvsactual_retained_n
    stats += [stats[15] / stats[0],                  #16: asvsactual_retained_p
              len(actualrejects),                    #17: asvsactual_rejected_n
              len(actualrejects) / stats[0]]         #18: asvsactual_rejected_p

    # Calculate score
    score = calc_score(stats[8], stats[1], stats[11], stats[2],
                       scoretype, weight)

    # Calculate estimates of 0: total input true targets, 1: total input true
    # nontargets, 2: total output true targets, 3: total output true nontargets
    estimates = estimate_true_values(stats[0], stats[3], stats[7], stats[1],
                                     stats[11], stats[2])

    # Store list of rejects or retains, whichever is shorter
    store = []
    if len(rejects) < 0.5 * len(asvs):
        store = ('reject', actualrejects)
    else:
        store = ('retain', asvs - actualrejects )

    # score, stats, estimates, hash, store
    rejects = tuple(sorted(actualrejects))
    return([score] + stats + estimates + [hash(rejects),  store])

def write_stats_and_cache(terms, thresholds, scores, filename, outdir):

    sh = open(os.path.join(outdir,f"{filename}_scores.csv"), 'w')
    ch = open(os.path.join(outdir,f"{filename}_resultcache"), 'w')

    # Write header
    head = ("score asvs_total targets_total nontargets_total "
            "asvsprelim_retained_n asvsprelim_retained_p "
            "asvsprelim_rejected_n asvsprelim_rejected_p "
            "targets_retained_n targets_retained_p "
            "targets_rejected_n targets_rejected_p "
            "nontargets_retained_n nontargets_retained_p "
            "nontargets_rejected_n nontargets_rejected_p "
            "asvsactual_retained_n asvsactual_retained_p "
            "asvsactual_rejected_n asvsactual_rejected_p "
            "inputtargets_total_estimate inputnontargets_total_estimate "
            "outputtargets_total_estimate "
            "outputtargets_total_estimate rejects_hash").split(" ")

    sh.write(",".join(["resultset"]
                      + [s + "_threshold" for s in terms]
                      + head) + '\n')


    # Write lines
    for i, (thresh, score) in enumerate(zip(thresholds, scores)):
        sh.write(",".join(str(v) for v in [i] + list(thresh) + score[:-1])
                 +'\n')
        ch.write("\t".join([str(i), score[-1][0]] + list(score[-1][1]))
                 +'\n')
    
    sh.close()
    ch.close()

def get_reject_from_store(asvs, n, store):
    if store[0] == 'retain':
        return(asvs - store[1])
    elif store[0] == 'reject':
        return(store[1])
    else:
        sys.exit(f"store type \'{store[0]}\' for resultset {n} is not "
                  "\'retain\' or \'reject\'")

def find_best_score(scores):
    #scores = stats
    
    scorelist = [ s[0] for s in scores ]
    sys.stdout.write(f"Minimum score of {min(scorelist)} achieved by "
                     f"{scorelist.count(min(scorelist))} threshold sets\n")

    return(sorted(scorelist))

def generate_resultsets(genval, stats, scoresort):
    #genval = args.generateASVresults
    i = 0
    if type(genval) is float:
        i = round(genval * len(scoresort)) - 1
    maxscore = scoresort[i]
    return([i for i, s in enumerate(stats) if s[0] <= maxscore])

def write_retained_asvs(infile, outfile, rejects):

    with open(infile, 'r') as infa, open(outfile, 'w') as outfa:
        for head, seq in SimpleFastaParser(infa):
            if head not in rejects:
                outfa.write(f">{head}\n{seq}\n")

def write_resultset_asvs(asvs, filename, infile, outdir, resultsets, store,
                         action, name):
    storei = 1 if action == 'dump' else -1
    
    for rs in resultsets:
        #rs = resultsets[0]
        out = f"{filename}_numtdumpresultset{rs}.fa" if name is None else name
        outfile = os.path.join(outdir, out)
        rsstore = store[rs][storei]
        rejects = get_reject_from_store(asvs, rs, rsstore)
        write_retained_asvs(infile, outfile, rejects)

def parse_resultcache(path, asvs):
    #path, asvs = [args.resultcache, set(raw['asvs'].keys())]
    fh = open(path, 'r')
    store = []
    for i, line in enumerate(fh):
        #i, line = next(enumerate(fh))
        vals = line.strip().split('\t')
        err = f"Error: {path} line {i+1}"
        if len(vals) < 2:
            sys.exit(f"{err} does not have enough items")
        setn = vals.pop(0)
        try:
            setn = int(setn)
        except:
            sys.exit(f"{err} starting value \'{setn}\' is not an integer")
        if setn != i:
                sys.exit(f"{err} starts with \'{setn}\', not \'{i}\'")
        action = vals.pop(0)
        if action not in {'reject', 'retain'}:
            sys.exit(f"{err} second value \'{action}\' is not 'reject' or "
                      "'retain'")
        missing = [a for a in vals if a not in asvs]
        if len(missing) > 0:
            sys.exit(f"{err} ASVs {', '.join(missing)} are not present or "
                      "identifiable in the supplied ASV file")
        store.append([setn, (action, vals)])
    fh.close()
    
    return(store)


