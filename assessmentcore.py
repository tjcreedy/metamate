#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Core functions for designating numts and optimising parameters"""

# Imports

import categorycounting

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

def calc_score(incorrect_known_good, total_known_good, incorrect_known_bad, total_known_bad, score_type, weight = 0.5):
	if(score_type == "standardised"):
		return(2 * weight * (incorrect_known_good / total_known_good) + 2 * (1-weight) * (incorrect_known_bad / total_known_bad))
	elif(score_type == "unstandardised"):
		return( ( weight * incorrect_known_good + (1 - weight) * incorrect_known_bad ) / ( weight * total_known_good + (1 - weight) * total_known_bad ) )
	else:
		sys.ext("Error: unknown score_type value passed to calc_score")

def check_good(failures, good):
	incorrect = failures.intersection(good)
	return(len(incorrect))

def check_bad(failures, bad):
	incorrect = bad - failures
	return(len(incorrect))

def resolve_ranges(range_string):
	# Split into segments
	ranges = re.split(",", range_string)
	# Separate parts of any range specifications
	ranges = [re.split("[-/]", r) for r in ranges]
	# TODO: add check for range specifications being complete
	# Convert any range specifications to full list
	ranges = [numpy.linspace(float(r[0]), float(r[1]), int(r[2])) if len(r) > 1 else [float(r[0]),] for r in ranges ]
	# Collapse list of lists
	return(list(itertools.chain(*ranges)))

def parse_spec(specfile):
	"Parse through the specification to expand terms and thresholds"
	
	# Open specifications file and remove comments
	
	with open(specfile) as fh:
		lines = fh.read().splitlines()
	
	commentline_filter = re.compile(r'^\s*$|^\s*#')
	lines_use = [l for l in lines if not commentline_filter.search(l)]
	lines_clean = [re.sub(r'#.*$', '', l) for l in lines_use]
	
	# Parse lines into specifications
	
	spec = {}
	thresh_list = [None] * len(lines_clean)
	thresh_total = 1
	spectext = [None] * len(lines_clean)
	
	for i, v in enumerate(lines_clean):
		
		# Clean and split up
		v = v.replace(" ","")
		values = re.split("\t+", v)
		
		# TODO: check there are 3 values!
		
		# Split terms into parts
		spectext[i] = values[0]
		values[0]  = re.split("[\|\+]", values[0])
		
		# TODO: check metric is n or p
		
		# Expand thresholds to list
		values[2] = resolve_ranges(values[2])
		
		# Add to list of thresholds
		thresh_total *= len(values[2])
		thresh_list[i] = values[2]
		
		# Add to dict
		spec[i] = dict(zip(['terms', 'metric', 'thresholds'], values))
	
	# TODO: Check thresh_total against maximum
	
	
	thresh_combos = list(itertools.product(*thresh_list))
	
	return(spectext, spec, thresh_combos)

def counts_from_spec(spec, data):
	"Does category counting for each individual specification in a specification set and returns a list of the category counts for each specification"
	
	# Set up new empty list for outputs
	counts = [None] * len(spec)
	
	# Work through specifications
	for specn, specdict in spec.items():
		
		# Check if partitioned and do counting
		if( len(specdict['terms']) == 1 ):
			counts[specn] = categorycounting.count_categories(data[specdict['terms'][0]], specdict['metric'] )
		else:
			partitiontermsdata = tuple( data[term] for term in specdict['terms'][1:] )
			counts[specn] = categorycounting.count_categories( categorycounting.multicategory( data[specdict['terms'][0]], partitiontermsdata), specdict['metric'] )
	
	return(counts)

def calc_combo_score(counts, thresholds, good, bad, score_type, weight = 0.5):
	"Finds failures and calculates score for a given set of category counts and a given set of thresholds. Counts and thresholds should be in two lists of equal lengths, where the nth item of the thresholds should be the threshold count for the nth counts"
	
	# Find all failures across counts
	failures = set(itertools.chain.from_iterable([ categorycounting.assess_failures(counts[i], t) for i, t in enumerate(thresholds) ]))
	
	# Get score
	score = calc_score(check_good(failures, good), len(good), check_bad(failures, bad), len(bad), score_type, weight)
	
	return(score, len(failures), len(good), len(bad), check_good(failures, good)/len(good), check_bad(failures, bad)/len(bad))

def get_minimum_thresholds(scores, threshold_combinations, spec):
	
	# Get list of scores only 
	score_list = [ score for score, failure_n, good_n, bad_n, good_fail_rate, bad_fail_rate in scores ]
	
	# Find minimum and indices
	minscore = min(score_list)
	indices = [i for i, v in enumerate(score_list) if v == minscore]
	
	# Generate thresholds and output text
	
	min_thresholds = []
	outtext = "Minimum score of " + str(minscore) + " achieved at:\n"
	n = 1
	
		# Go through each combination of thresholds resulting in an minimum
	for min_index in indices:
		
		threshold_comb = threshold_combinations[min_index]
		min_thresholds.append(threshold_comb)
		
		outtext += "\t(" + str(n) + ")\n"
		
		# Go through each threshold
		for i, v in enumerate(threshold_comb):
			outtext += "\tThreshold "
			
			if(spec[i]['metric'] == "n"):
				outtext += "number "
			else:
				outtext += "proportion "
			outtext += "of " + spec[i]['terms'][0] + " read counts per haplotype"
			
			if(len(spec[i]['terms']) > 1):
				outtext += " by " + " and ".join(spec[i]['terms'][1:])
			
			outtext += " = " + str(v) + "\n"
		
		outtext += "\tPreliminary number of putative exclusions: " + str(scores[min_index][1]) + "\n"
		outtext += "\tProportion of \'good\' haplotypes incorrectly excluded: " + str(scores[min_index][4]) + "\n"
		outtext += "\tProportion of \'bad\' haplotypes incorrectly included: " + str(scores[min_index][5]) + "\n"
		outtext += "\tFinal number of exclusions: " + str(int(scores[min_index][1] + scores[min_index][5] * scores[min_index][3] - scores[min_index][4] * scores[min_index][2])) + "\n"
		
		n += 1
	
	return(minscore, min_thresholds, outtext)

def output_filtered_haplotypes(counts, min_thresholds, good, bad, file, filename, outdir):
	
	#if(type(min_thresholds[0]) != list):
	#	min_thresholds = [min_thresholds,]
	
	#min_i, thresholds = list(enumerate(min_thresholds))[10]
	
	for min_i, thresholds in enumerate(min_thresholds):
		
		failures = set(itertools.chain.from_iterable([ categorycounting.assess_failures(counts[i], t) for i, t in enumerate(thresholds) ]))
		
		exclude = failures.union(bad) - good
		
		with open(os.path.join(outdir, filename + "_filtered_set" + str(min_i) + ".fa"), "w") as outfasta:
			
			for head, seq in SimpleFastaParser(file):
			
				if( head not in exclude):
					outfasta.write(">%s\n%s\n" % (head, seq))
