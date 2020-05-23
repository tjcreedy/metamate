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

def parse_specfile(args):
    "Parse through the specification to expand terms and thresholds"
    
    filename = os.path.splitext(os.path.basename(args.asvs))[0]
    
    # Open specifications file and parse lines into specifications
    fh = open(args.specification, 'r')
    spectext = ''
    for l in fh:
        #l = fh.readline()
        if re.match("^\s*#|^$", l):
            continue
        else:
            spectext += l.strip()
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
            thresholds.append(specs[term] if term in specs else [float('nan')])
        threshcombos.extend(itertools.product(*thresholds))
    # Check 
    if( not args.taxgroups and any("taxon" in e for e in spectext)):
        sys.exit("Error, taxon specified as a binning strategy but no taxon "
                 "file supplied")
    #Output specifications in table
    path = os.path.join(args.outputdirectory, f"{filename}_specifications.csv")
    with open(path, 'w') as o:
        o.write(",".join(specdict.keys()) + "\n")
        for thresh in threshcombos:
            o.write(",".join([str(t) for t in thresh]) + "\n")
    
    return(specdict, termorder, threshcombos)

def get_validated(raw, minmaxbp, args, filename):
    
    # Create length-based control list
    sys.stdout.write("Identifying control non-target ASVs based on length.\n")
    
    nontargetlength = filterlength.check_length_multi(raw['asvs'], minmaxbp,
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
    sys.stdout.write("Identifying control target ASVs based on references.\n")
    
    refmatch = filterreference.blast_for_presence(raw['path'],
                                                   args.outputdirectory,
                                                   args.reference,
                                                   args.matchpercent,
                                                   args.matchlength,
                                                   args.threads)
    
    target = refmatch - nontarget
    
    # Finalise targets
    if len(target) > 0:
        sys.stdout.write(f"Found {len(target)} target ASVs: {len(refmatch)} "
                          "out of all ASVs matched to reference set, "
                         f"{len(refmatch - target)} non-target ASVs removed."
                         "\n")
    else:
        err = "Error: no target ASVs found"
        if len(refmatch) > 0:
            err = (f"{err}, although {len(refmatch)} ASVs matched to "
                    "reference set. Length thresholds may be too stringent "
                    "and find too many non-target ASVs.")
        else:
            err = (f"{err}. Check your reference and your reference matching "
                    "thresholds.")
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

def calc_stats(counts, termorder, asvs, target, nontarget, anythreshold, 
               scoretype, thresholds, weight = 0.5):
    #asvs, anythreshold, scoretype, thresholds, weight =  [set(raw['asvs'].keys()), args.anythreshold, "standardised" , thresholdcombos[0], 0.5]
    """For a given set of category counts and a given set of thresholds, counts
    retention, calculates scores and estimates statistics. Counts and 
    thresholds should be in two lists of equal lengths, where the nth item of 
    the thresholds should be the threshold count for the nth counts"""
    
    # Find all rejected_asvs across counts
    rejects = []
    for term, thresh in zip(termorder, thresholds):
        #term, thresh = list(zip(termorder, thresholds))[3]
        rejects.extend(categorycounting.reject(counts[term], thresh,
                                               anythreshold))
    rejects = set(rejects)
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
    stats += [stats[13 / stats[2]],                  #14: nontargets_rejected_p
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
        ('rejects', actualrejects)
    else:
        ('retain', asvs - actualrejects )
    
    # score, stats, estimates, hash, store
    rejects = tuple(sorted(actualrejects))
    return([score] + stats + estimates + [hash(rejects),  store])

def write_specs_and_stats(terms, thresholds, scores, path):
    
    with open(path, "w") as o:
        
        # Write header
        head = ("score targets_total nontargets_total asvsprelim_retained_n "
                "asvsprelim_retained_p asvsprelim_rejected_n "
                "asvsprelim_rejected_p targets_retained_n targets_retained_p "
                "targets_rejected_n targets_rejected_p nontargets_retained_n "
                "nontargets_retained_p nontargets_rejected_n "
                "nontargets_rejected_p asvsactual_retained_n "
                "asvsactual_retained_p asvsactual_rejected_n "
                "asvsactual_rejected_p inputtargets_total_estimate "
                "inputnontargets_total_estimate outputtargets_total_estimate "
                "outputtargets_total_estimate").split(" ")
        
        o.write(",".join([s + "_threshold" for s in terms] + head) + '\n')
        
        # Write lines
        for thresh, score in zip(thresholds, scores):
            o.write(",".join(str(v) for v in list(thresh) + score[:-2])+'\n')

def get_minimum_thresholds(scores, terms, thresholdcombinations, spec):
    #scores, thresholdcombinations, spec = [stats, thresholdcombos, specs]
    # Get list of scores only
    scorelist = [ s[0] for s in scores ]
    
    # Find minimum and indices
    minscore = min(scorelist)
    indices = [i for i, v in enumerate(scorelist) if v == minscore]
    
    # Generate thresholds and output text
    
    minthresholds = []
    outtext = f"Minimum score of {str(minscore)} achieved at:\n"
    n = 1
    
    # Go through each combination of thresholds resulting in an minimum
    for minindex in indices:
        #minindex = indices[0]
        
        threshcomb = thresholdcombinations[minindex]
        minthresholds.append(threshcomb)
        
        outtext += f"\t({str(n)})\n"
        
        # Go through each threshold
        for term, thresh in zip(terms, threshcomb):
            #i, v = list(enumerate(threshcomb))[0]
            outtext += "\tThreshold "
            
            if(spec[term]['metric'] == "n"):
                outtext += "number "
            else:
                outtext += "proportion "
            outtext += f"of {spec[term]['terms'][0]} read counts per haplotype"
            
            if(len(spec[term]['terms']) > 1):
                outtext += " by " + " and ".join(spec[terms]['terms'][1:])
            
            outtext += f" = {str(thresh)}\n"
        
        outtext += ( "\tPreliminary number of putative exclusions: "
                    f"{str(scores[minindex][4])}\n\t"
                     "Proportion of \'verified target\' ASVs retained: "
                    f"{str(scores[minindex][6]/scores[minindex][2])}\n\t"
                     "Proportion of \'verified non-target\' ASVs retained: "
                    f"{str(scores[minindex][7]/scores[minindex][3])}\n\t"
                    f"Final number of exclusions: {str(scores[minindex][8])}\n"
                    )
        n += 1
    
    return(minscore, minthresholds, outtext)


def output_filtered_haplotypes(counts, minthresholds, termorder, anythreshold,
                               good, bad, file, filename, outdir):
    # minthresholds, anythreshold, good, bad, file, outdir = [minthresh, args.anythreshold, target, nontarget, raw['path'], args.outputdirectory]
    
    for mini, thresholds in enumerate(minthresholds):
        #mini, thresholds = list(enumerate(minthresholds))[1]
        rejects = []
        for term, thresh in zip(termorder, thresholds):
            rejects.extend(categorycounting.reject(counts[term], term, 
                                                   anythreshold))
        rejects = set(rejects)
        
        exclude = rejects.union(bad) - good
        with open(file) as infasta:
            filen = mini + 1
            outname = os.path.join(outdir, 
                                   f"{filename}_filtered_set{str(filen)}.fa") 
            with open(outname, "w") as outfasta:
                for head, seq in SimpleFastaParser(infasta):
                    if( head not in exclude):
                        outfasta.write(">%s\n%s\n" % (head, seq))
