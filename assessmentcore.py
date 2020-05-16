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

def parse_spec(args):
    "Parse through the specification to expand terms and thresholds"
    
    filename = os.path.splitext(os.path.basename(args.asvs))[0]
    
    # Set up specifications
    spec = {}
    threshlist = []
    threshtotal = 1
    spectext = []
    
    # Open specifications file and parse lines into specifications
    fh = open(args.specification, 'r')
    linen = 0
    specn = 0
    for l in fh.readline():
        l = fh.readline()
        l.strip()
        linen += 1
        if re.match("^\s*#|^$", l):
            continue
        else:
            # Clean and split up
            v = l.replace('#.*$| ', '')
            values = re.split("\t+", v)
            
            # Set up error
            err = (f"Error, malformed specification in {args.specification} "
                   f"line {linen}")
            
            # Check there are 3 values
            if len(values) != 3:
                sys.exit(f"{err}: specification line should have three "
                         "tab-separated entries")
            
            # Split terms into parts
            spectext.append(f"{values[0]}({values[1]})")
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
            # Add null value if requested
            if(args.addnull):
                if(values[1] == 'p' and not 1 in values[2]):
                    values[2].append(1)
                elif(values[1] == 'n'):
                    values[2].append(99999999)
            
            # Add to list of thresholds
            threshtotal *= len(values[2])
            threshlist.append(values[2])
            
            # Add to dict
            spec[specn] = dict(zip(['terms', 'metric', 'thresholds'], values))
            specn += 1
    
    # TODO: Check thresh_total against maximum
    threshcombos = list(itertools.product(*threshlist))
    
    # Check 
    if( not args.taxgroups and any("taxon" in e for e in spectext)):
        sys.exit("Error, taxon specified as a binning strategy but no taxon "
                 "file supplied")
    
    #Output specifications in table
    path = os.path.join(args.outputdirectory, f"{filename}_specifications.csv")
    with open(path, 'w') as o:
        o.write(",".join(spectext) + "\n")
        for thresh in threshcombos:
            o.write(",".join([str(t) for t in thresh]) + "\n")
    
    return(spectext, spec, threshcombos)

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
                         f"{len(refmatch-nontarget)} non-target ASVs removed."
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
        for asv, cat in zip(["lengthfail", "stopfail", "refpass"], 
                            [nontargetlength, nontargettrans, refmatch]):
            for a in asv: o.write(f"{cat}\t{a}\n")
    
    return(target, nontarget)


def counts_from_spec(spec, data):
    """Count categories for each individual specification in a specification 
    set and return a list of the category counts for each specification"""
    # spec = specs
    
    
    counts = []
    
    # Work through specifications
    for specn, specdict in spec.items():
        #specn, specdict = list(spec.items())[1]
        # Check if partitioned and do counting
        if len(specdict['terms']) == 1 :
            counts.append(categorycounting.count_categories(
                                                    data[specdict['terms'][0]],
                                                    specdict['metric']))
        else:
            parttermsdata = tuple( data[t] for t in specdict['terms'][1:] )
            multicat = categorycounting.multicategory(
                                                    data[specdict['terms'][0]],
                                                    parttermsdata)
            counts.append(categorycounting.count_categories(multicat,
                                                           specdict['metric']))
    
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

def calc_stats(counts, asvs, target, nontarget, anythreshold, 
               scoretype, thresholds, weight = 0.5):
    #asvs, anythreshold, score_type, thresholds, weight =  [set(raw['asvs'].keys()), args.anythreshold, "standardised" , thresholdcombos[467], 0.5]
    """For a given set of category counts and a given set of thresholds, counts
    retention, calculates scores and estimates statistics. Counts and 
    thresholds should be in two lists of equal lengths, where the nth item of 
    the thresholds should be the threshold count for the nth counts"""
    
    # Find all rejected_asvs across counts
    rejects = []
    for i, t in enumerate(thresholds):
        #i, t = list(enumerate(thresholds))[1]
        rejects.extend(categorycounting.reject(counts[i], t, anythreshold))
    rejects = set(rejects)
    
    # Input counts
    inputs = [len(asvs), len(target), len(nontarget), len(rejects)]
    inputs.append(inputs[0]-inputs[3])
    
    # Calculate number of retained target and nontarget asvs, plus number 
    # of actual retentions 
    retained_vals = [len(target - rejects), len(nontarget - rejects),
                     len(rejects - target)]
    
    # Calculate score
    score = calc_score(retained_vals[0], inputs[1], retained_vals[1], 
                       inputs[2], scoretype, weight)
    
    # Calculate estimates of 0: total input true targets, 1: total input true
    # nontargets, 2: total output true targets, 3: total output true nontargets
    estimates = estimate_true_values(inputs[0], inputs[4], retained_vals[0],
                                     inputs[1], retained_vals[1], inputs[3])
    
    # score, asvs, target, nontarget, rejectedasvs, retainedasvs, 
    # retained_target, retained_nontarget, actual_retainedasvs, true_target, 
    # true_nontarget, true_retained_target, true_retained_nontarget,
    # rejectedasvshash
    rejects = tuple(sorted(rejects))
    return([score] + inputs + retained_vals + estimates + [hash(rejects)])

def write_specs_and_stats(specs, thresholds, scores, path):
    
    with open(path, "w") as o:
        
        # Write header
        scorehead = ["score", "asvs", "target", "nontarget", "rejectedasvs",
                     "retainedasvs", "retained_target", "retained_nontarget",
                     "actual_retainedasvs", "est_true_target",
                     "est_true_nontarget", "est_true_retained_target",
                     "est_true_retained_nontarget, hash_rejectedasvs"]
        
        o.write(",".join([s + "_threshold" for s in specs] + scorehead) + '\n')
        
        # Write lines
        for thresh, score in zip(thresholds, scores):
            o.write(",".join(str(v) for v in list(thresh) + score)+'\n')


def get_minimum_thresholds(scores, thresholdcombinations, spec):
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
        for i, v in enumerate(threshcomb):
            #i, v = list(enumerate(threshcomb))[0]
            outtext += "\tThreshold "
            
            if(spec[i]['metric'] == "n"):
                outtext += "number "
            else:
                outtext += "proportion "
            outtext += f"of {spec[i]['terms'][0]} read counts per haplotype"
            
            if(len(spec[i]['terms']) > 1):
                outtext += " by " + " and ".join(spec[i]['terms'][1:])
            
            outtext += f" = {str(v)}\n"
        
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


def output_filtered_haplotypes(counts, minthresholds, anythreshold, good, bad,
                               file, filename, outdir):
    
    #if(type(min_thresholds[0]) != list):
    #    min_thresholds = [min_thresholds,]
    
    #min_i, thresholds = list(enumerate(min_thresholds))[10]
    
    for mini, thresholds in enumerate(minthresholds):
        rejectedasvs = set()
        for i, t in enumerate(thresholds):
            rejectedasvs.add(categorycounting.reject(counts[i], t, 
                                                     anythreshold))
        
        exclude = rejectedasvs.union(bad) - good
        with open(file) as infasta:
            filen = mini + 1
            filename = os.path.join(outdir, 
                                   f"{filename}_filtered_set{str(filen)}.fa") 
            with open(filename, "w") as outfasta:
                for head, seq in SimpleFastaParser(infasta):
                    if( head not in exclude):
                        outfasta.write(">%s\n%s\n" % (head, seq))
