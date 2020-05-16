#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Functions for assessing read counts by category combinations against a threshold"""

# Imports

from collections import defaultdict
import sys

# Global variables

# Class definitions

# Function definitions

# Categories will be stored as a dict of the form

    #dict[category][name] : count    for reads by library
    #dict[name] : count        for read totals by zotu
    #dict[category] : name        for clades and taxa

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

def reject(counts, threshold, anythreshold):
    #counts, threshold = counts[i], t
    # Work through sets of scores
    #name = "uniq13;size=257"
    #count = counts[name]
    out = []
    for name, count in counts.items():
        if ((anythreshold and min(count) < threshold)
            or (not anythreshold and all(c < threshold for c in count))):
            out.append(name)
    return(out)

def multicategory(countdict, otherdicts):
    
    """Find all combinations of sets in category dictionaries.
    First dictionary must have counts of incidences, i.e {x : {a : 1, b : 2}},
    the rest are passed in a tuple of length >= 1, each can have counts or 
    just sets, i.e. {y : {a , b}}
    """
    
    # Check otherdicts is a tuple or error
    
    if(type(otherdicts) != tuple):
        otherdicts = (otherdicts,)
    
    # Initialise master dict with count dict
    multidict = countdict
    
    # Work through each further dict combining with master
    for newdict in otherdicts:
        
        # Create new multdict
        newmulti = dict()
        
        # Work through categories in master
        for mcat, mcounts in multidict.items():
            
            # Work through categories in new dict
            for ncat, nvalue in newdict.items():
                
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



#x = {"l1" : { "A" : 3, "B" : 2}, "l2" : { "B" : 1, "C" : 3}, "l3" : { "A" : 2, "C" : 1, "D" : 4}}
#y = {"t1" : { "A" : 5}, "t2" : { "B" : 3, "C" : 4}, "t3" : {"D" : 4}}
#z = {"c1" : { "A" , "B" }, "c2" : {"C"}, "c3" : {"D"}}
#xyz = multicategory(x,(y,z))
#path = "test.csv"
#asvs = ["A", "B", "C", "D"]
#count_dict = x
