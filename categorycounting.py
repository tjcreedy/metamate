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

	#dict[category][name] : count	for reads by library
	#dict[name] : count		for read totals by zotu
	#dict[category] : name		for clades and taxa

def count_categories(catdict, metric):
	#catdict, metric = [data[specdict['terms'][0]], specdict['metric']]
	
	# Set up dictiony of sets for scores
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
	reject = {name for name, count in counts.items() if (anythreshold and min(count) < threshold) or (not anythreshold and all(c < threshold for c in count))}
	
	return(reject)


def multicategory(count_dict, other_dicts_tuple):
	
	"Finds all combinations of sets in category dictionaries. First dictionary must have counts of incidences, i.e {x : {a : 1, b : 2}}, rest are passed in a tuple of length >=, each can have counts or just sets, i.e. {y : {a , b}}"
	
	# Check other_dicts_tuple is a tuple or error
	
	if(type(other_dicts_tuple) != tuple):
		other_dicts_tuple = (other_dicts_tuple,)
	
	# Initialise master dict with count dict
	multi_dict = count_dict
	
	# Work through each further dict combining with master
	for add_dict in other_dicts_tuple:
		
		# Create new multdict
		new_multi = dict()
		
		# Work through categories in master
		for mcat, mcounts in multi_dict.items():
			
			# Work through categories in new dict
			for acat, avalue in add_dict.items():
				
				# Check if avalue is set of names or dictionary of name:counts
				avalue_set = avalue
				if type(avalue) == "dict":
					avalue_set = avalue.keys()
				
				# Find all names shared between current categories
				nnames = set( mcounts.keys() ).intersection( avalue_set )
				
				# Skip combination if no names shared
				
				if len(nnames) == 0:
					continue
				
				# Create the combined category name
				ncat = mcat+acat
				
				# Create a combined subdictionary
				new_multi[ncat] = dict()
				
				# Add the minimum value of the two counts for each shared name
				if type(avalue) == "dict":
					new_multi[ncat] = {nname : min(mcounts[nname], avalue[nname]) for nname in nnames}
				else:
					new_multi[ncat] = {nname : mcounts[nname] for nname in nnames}
				
			
		
		# Overwrite old multi_dict
		multi_dict = new_multi
	
	return(multi_dict)

def write_count_dict(count_dict, asvs, path):
	""
	
	# TODO: Check the count dict is of the right structure
	
	# Order asvs name
	
	asvs_sort = sorted(asvs)
	
	# Write to file
	
	with open(path, 'w') as o:
		
		# Write header
		
		o.write(",".join([""]+asvs_sort)+"\n")
		
		# Write data lines
		
		for lib, counts in count_dict.items():
			values = [str(count_dict[lib][asv]) if asv in count_dict[lib].keys() else str(0) for asv in asvs_sort]
			o.write(",".join([lib]+values)+"\n")



#x = {"l1" : { "A" : 3, "B" : 2}, "l2" : { "B" : 1, "C" : 3}, "l3" : { "A" : 2, "C" : 1, "D" : 4}}
#y = {"t1" : { "A" : 5}, "t2" : { "B" : 3, "C" : 4}, "t3" : {"D" : 4}}
#z = {"c1" : { "A" , "B" }, "c2" : {"C"}, "c3" : {"D"}}
#xyz = multicategory(x,(y,z))
#path = "test.csv"
#asvs = ["A", "B", "C", "D"]
#count_dict = x