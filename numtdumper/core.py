#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core functions for designating numts and optimising parameters"""

# Imports

import os
import sys
import shutil
import re
import functools
import gzip
import numpy
import itertools

from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict

from numtdumper import filterlength
from numtdumper import filtertranslate
from numtdumper import filterreference

def count_categories(catdict, metric):
    #catdict, metric = [data[specdict['terms'][0]], specdict['metric']]
    
    # Set up dict of sets for scores
    counts = defaultdict(set)
    
    # Work through each category
    for category, catcounts in catdict.items():
        #category, catcounts = list(catdict.items())[0]
        # Get total for category
        total = sum(catcounts.values())
        
        # Work through each read
        for name, count in catcounts.items():
            # Get proportion for this read and add to dict
            if(metric == "p"):
                counts[name].add(count/total)
            elif(metric == "n"):
                counts[name].add(count)
            else:
                sys.exit("Error: unknown metric passed to score_categories")
        
    
    return(counts)

def reject(counts, threshold, anyfail):
    #counts, threshold = counts[i], t
    # Work through sets of scores
    #name = "uniq13;size=257"
    #count = counts[name]
    out = []
    for name, count in counts.items():
        if ((anyfail and min(count) < threshold)
            or (not anyfail and all(c < threshold for c in count))):
            out.append(name)
    return(out)

def multicategory(countdict, otherdicts):
    
    """Find all combinations of sets in category dictionaries.
    First dictionary must have counts of incidences, i.e {x : {a : 1, b : 2}},
    the rest are passed in a tuple of length >= 1, each can have counts or 
    just sets, i.e. {y : {a , b}}
    """
    # countdict, otherdicts = [data[specdict['terms'][0]], parttermsdata]
    # Check otherdicts is a tuple or error
    
    if(type(otherdicts) != tuple):
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
                    newmulti[ncat] = {nn : min(mcounts[nn], nvalue[nn]) 
                                      for nn in nnames}
                else:
                    newmulti[ncat] = {nn : mcounts[nn] for nn in nnames}
                
            
        
        # Overwrite old multi_dict
        multidict = newmulti
    
    return(multidict)

def write_count_dict(countdict, asvs, path):
    ""
    
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
            rangeout.append(numpy.linspace(float(r[0]),
                                           float(r[1]),
                                           int(r[2])))
        else:
            rangeout.append([float(r[0])])
    # Collapse list of lists
    return(list(itertools.chain(*rangeout)))

def resolve_spec(spec):
    
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

def parse_specs(args, null = float('nan')):
    "Parse through the specification to expand terms and thresholds"
    #null = 0
    
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
    
    # Check
    if( not args.taxgroups and any("taxon" in e for e in spectext)):
        sys.exit("Error, taxon specified as a binning strategy but no taxon "
                 "file supplied")
    
    # Clean up and parse for additive and multiplicative specs
    spectext = re.sub(' ', '', spectext)
    specslist = [s.split('*') for s in spectext.split('+')]
    
    # Find unique terms and extract thresholds
    specdict = dict()
    termvals = []
    for specl in specslist:
        out = dict()
        for spec in specl:
            l = resolve_spec(spec)
            specdict[l[0]] = l[1:3]
            out[l[0]] = l[3]
        termvals.append(out)
    
    specs = {'name': list(specdict.keys())}
    specs.update({k: v for k, v in zip(['spec', 'metric'], 
                                       map(list, zip(*specdict.values())))})
    
    # Build all threshold combinations
    threshlists = []
    terms = []
    nthresh = 0
    for tval in termvals:
        #tval = termvals[6]
        specl = []
        threshl = []
        for spec in specs['name']:
            threshl.append(tval[spec] if spec in tval else [null])
            if spec in tval: specl.append(spec)
        n = functools.reduce(lambda x, y: x * y, 
                                   [len(t) for t in threshl], 1)
        nthresh += n
        terms.append([specl, nthresh])
        threshlists.append(threshl)
    nterm = len(terms)
    
    if args.action == 'dump' and nthresh > 1:
        sys.exit("Error: mode 'dump' allows only a single threshold set. To "
                 "run wih multiple threshold sets, use mode 'find'\n")
    
    # Construct generators
    terms = (l for t, c in terms for l in [t] * c)
    thresholds = (t for tl in threshlists for t in itertools.product(*tl))
    
    return(specs, nterm, nthresh, terms, thresholds)

def get_validated(raw, args, filename):
    #filname = infilename
    
    # Create length-based control list
    sys.stdout.write("Identifying control non-target ASVs based on length...")
    
    nontargetlength = filterlength.check_length_multi(raw['asvs'],
                                                      [args.minimumlength,
                                                      args.maximumlength],
                                                      args, fail=True)
    
    sys.stdout.write(f"found {len(nontargetlength)} candidates\n")
    
    # Create translation based control list
    sys.stdout.write("Identifying control non-target ASVs based on "
                     "translation...")
    
    nontargettrans = filtertranslate.check_stops_multi(raw['asvs'], args,
                                                       fail = True)
    
    sys.stdout.write(f"found {len(nontargettrans)} candidates\n")
    
    # Finalise nontargets
    nontarget = set(nontargetlength + nontargettrans)
    overlap = set(nontargetlength).intersection(set(nontargettrans))
    if len(nontarget) > 0:
        sys.stdout.write(f"Designating a total of {len(nontarget)} unique "
                         f"ASVs as non-target ({len(overlap)} were found in "
                          "both sets)\n")
    else:
        sys.exit("Error: no non-target ASVs could be found. Prior data "
                 "filtering may have been too stringent.")
    
    # Create reference based control list
    
    wd = filterreference.make_temp_blastwd(args.outputdirectory, "blastdb")
    allcandidates = []
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
        allcandidates.extend(candidates)
    
    shutil.rmtree(wd)
    
    refmatch = set(allcandidates)
    target = refmatch - nontarget
    
    # Finalise targets
    if len(target) > 0:
        sys.stdout.write(f"Designating a total of {len(target)} unique ASVs "
                         f"as target ({len(refmatch)} unique matched to "
                          "references and/or blast database, "
                         f"{len(refmatch - target)} rejected due to inclusion "
                          "in non-target set)\n")
    else:
        err = "Error: no target ASVs found"
        if len(refmatch) > 0:
            err = (f"{err}, although {len(refmatch)} ASVs matched to "
                    "reference set. Length thresholds may be too stringent "
                    "and find too many non-target ASVs, or reference dataset "
                    "is not sufficiently curated")
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

def counts_from_spec(specs, data):
    """Count categories for each individual specification in a specification
    set and return a list of the category counts for each specification"""
    
    counts = dict()
    
    # Work through specifications
    for name, spec, metric in zip(*specs.values()):
        #name, spec, metric = list(zip(*specs.values()))[1]
        # Check if partitioned and do counting
        if len(spec) == 1 :
            counts[name] = count_categories(data[spec[0]],
                                                             metric)
        else:
            partspecdata = tuple( data[s] for s in spec[1:] )
            multicat = multicategory(data[spec[0]],
                                                      partspecdata)
            counts[name] = count_categories(multicat, metric)
    
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

def assess_numts(threshnames, counts, anyfail,
                 asvs, target, nontarget, scoretype, weight,
                 thresholds):
    rejects = apply_reject(threshnames, thresholds, counts, anyfail)
    stats = calc_stats(rejects, asvs, target, nontarget, scoretype, weight)
    
    # TODO: send thresholds + stats to a parallel writer instead, and only
    # output sufficient information to generate resultsets
    
    return(list(thresholds) + stats)

def apply_reject(threshnames, thresholds, counts, anyfail):
    rejects = []
    for term, thresh in zip(threshnames, thresholds):
        #term, thresh = list(zip(threshnames, thresholds))[3]
        rejects.extend(reject(counts[term], thresh, anyfail))
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

def write_stats_and_cache(specs, scores, terms, filename, outdir):
    
    sh = open(os.path.join(outdir,f"{filename}_results.csv"), 'w')
    ch = gzip.open(os.path.join(outdir,f"{filename}_resultcache"), 'wt')
    
    # Write header
    head = ("score asvs_total targets_total_observed nontargets_total_observed"
            " asvsprelim_retained_n asvsprelim_retained_p "
            "asvsprelim_rejected_n asvsprelim_rejected_p "
            "targets_retained_n targets_retained_p "
            "targets_rejected_n targets_rejected_p "
            "nontargets_retained_n nontargets_retained_p "
            "nontargets_rejected_n nontargets_rejected_p "
            "asvsactual_retained_n asvsactual_retained_p "
            "asvsactual_rejected_n asvsactual_rejected_p "
            "targets_total_estimate nontargets_total_estimate "
            "targets_retained_estimate nontargets_retained_estimate "
            "rejects_hash").split(" ")
    
    sh.write(",".join(["resultindex", "term"]
                      + [s + "_threshold" for s in specs['name']]
                      + head) + '\n')
    
    # Write lines
    for i, term, score in zip(itertools.count(), terms, scores):
        sh.write(",".join(str(v) for v in [i, '*'.join(term)] + score[:-1])
                 +'\n')
        ch.write('\t'.join([str(i), score[-1][0]] + list(score[-1][1])) +'\n')
    
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

def find_best_score(scores, nspec):
    #scores, nspec = stats, len(specs['name'])
    
    scorelist = [ s[nspec] for s in scores ]
    sys.stdout.write(f"Minimum score of {min(scorelist)} achieved by "
                     f"{scorelist.count(min(scorelist))} threshold sets\n")
    
    return(sorted(scorelist))

def generate_resultsets(genval, stats, scoresort, nspec):
    #genval = args.generateASVresults
    i = 0
    if type(genval) is float:
        i = round(genval * len(scoresort)) - 1
    maxscore = scoresort[i]
    return([i for i, s in enumerate(stats) if s[nspec] <= maxscore])

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
    fh = gzip.open(path, 'rt')
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


