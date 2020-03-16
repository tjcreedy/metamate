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

def parse_spec(specfile, addnull):
	"Parse through the specification to expand terms and thresholds"
	
	# Open specifications file and remove comments
	
	with open(specfile) as fh:
		lines = fh.read().splitlines()
	
	commentline_filter = re.compile(r'^\s*$|^\s*#')
	
	# Parse lines into specifications
	spec = {}
	thresh_list = []
	thresh_total = 1
	spectext = []
	i = 0
	
	for l in lines:
		if(commentline_filter.search(l)):
			next
		else:
			# Clean and split up
			v = l.replace('#.*$| ', '')
			values = re.split("\t+", v)
			
			# Set up error
			err = "Error, malformed specification in " + specfile + " line " + str(i+1)
			
			# Check there are 3 values
			if(len(values) != 3):
				errlen = err + ": specification line should have three tab-separated entries"
				sys.exit(errlen)
			
			# Split terms into parts
			spectext.append(values[0]+"("+values[1]+")")
			values[0]  = re.split("[\|\+]", values[0])
			
			# Check terms are formed correctly
			if(not values[0][0] in ['library', 'total']):
				errterm = err + ": first term must be \'library\' or \'total\'"
				sys.exit(errterm)
			
			if(len(values[0]) > 1 and not all(t in ['taxon', 'clade'] for t in values[0][1:])):
				errterm = err + ": secondary term(s) must be \'clade\' or \'taxon\'"
				sys.exit(errterm)
			
			# Check metric is n or p
			if(not values[1] in ['n', 'p']):
				errmet = err = ": metric should be \'n\' or \'p\'"
				sys.exit(errmet)
			
			# Expand thresholds to list
			values[2] = resolve_ranges(values[2])
			
			# Add null value if requested
			if(addnull):
				if(values[1] == 'p' and not 1 in values[2]):
					values[2].append(1)
				elif(values[1] == 'n'):
					values[2].append(99999999)
			
			# Add to list of thresholds
			thresh_total *= len(values[2])
			thresh_list.append(values[2])
			
			# Add to dict
			spec[i] = dict(zip(['terms', 'metric', 'thresholds'], values))
			i += 1
	
	# TODO: Check thresh_total against maximum
	
	thresh_combos = list(itertools.product(*thresh_list))
	
	return(spectext, spec, thresh_combos)

def counts_from_spec(spec, data):
	# spec = specs
	"Does category counting for each individual specification in a specification set and returns a list of the category counts for each specification"
	
	# Set up new empty list for outputs
	counts = [None] * len(spec)
	
	# Work through specifications
	for specn, specdict in spec.items():
		#specn, specdict = list(spec.items())[1]
		# Check if partitioned and do counting
		if( len(specdict['terms']) == 1 ):
			counts[specn] = categorycounting.count_categories(data[specdict['terms'][0]], specdict['metric'])
		else:
			partitiontermsdata = tuple( data[term] for term in specdict['terms'][1:] )
			counts[specn] = categorycounting.count_categories( categorycounting.multicategory( data[specdict['terms'][0]], partitiontermsdata), specdict['metric'] )
	
	return(counts)

def calc_score(retained_target, target, retained_nontarget, nontarget, score_type, weight = 0.5):
	if(score_type == "standardised"):
		return(2 * ( weight * (1 - retained_target/target) + (1-weight) * retained_nontarget/nontarget))
	elif(score_type == "unstandardised"):
		return( ( weight * (target - retained_target) + (1 - weight) * (nontarget - retained_target)) / ( weight * target + (1 - weight) * nontarget) )
	else:
		sys.ext("Error: unknown score_type value passed to calc_score")

def estimate_true_values(asvs, retained_asvs, retained_target, target, retained_nontarget, nontarget):
	try:
		true_target = ( retained_asvs - asvs * (retained_nontarget / nontarget) ) / (retained_target / target - retained_nontarget / nontarget )
		true_nontarget = asvs - true_target
		true_retained_target = true_target * (retained_target / target)
		true_retained_nontarget = (retained_nontarget / nontarget) * (asvs - true_target)
		out = [true_target, true_nontarget, true_retained_target, true_retained_nontarget]
	except:
		out = ['NA', 'NA', 'NA', 'NA']
	return(out)

def calc_stats(counts, asvs, target, nontarget, anythreshold, score_type, thresholds, weight = 0.5):
	#asvs, anythreshold, score_type, thresholds, weight =  [set(rawasvs.keys()), args.anythreshold, "standardised" , threshold_combinations[467], 0.5]
	"For a given set of category counts and a given set of thresholds, counts retention, calculates scores and estimates statistics. Counts and thresholds should be in two lists of equal lengths, where the nth item of the thresholds should be the threshold count for the nth counts"
	
	# Find all rejected_asvs across counts
	#i, t = list(enumerate(thresholds))[1]
	rejectedasvs = set(itertools.chain.from_iterable([ categorycounting.reject(counts[i], t, anythreshold) for i, t in enumerate(thresholds) ]))
	
	# Input counts
	inputs = [len(asvs), len(target), len(nontarget), len(rejectedasvs)]
	inputs.append(inputs[0]-inputs[3])
	
	# Calculate number of retained target and nontarget asvs, plus number of actual retentions 
	retained_vals = [len(target - rejectedasvs), len(nontarget - rejectedasvs), len(rejectedasvs - target)]
	
	# Calculate score
	score = calc_score(retained_vals[0], inputs[1], retained_vals[1], inputs[2], score_type, weight)
	
	# Calculate estimates of total input true targets, total input true nontargets, total output true targets, total output true nontargets
	estimates = estimate_true_values(inputs[0], inputs[4], retained_vals[0], inputs[1], retained_vals[1], inputs[3])
	
	# score, asvs, target, nontarget, rejectedasvs, retainedasvs, retained_target, retained_nontarget, actual_retainedasvs, true_target, true_nontarget, true_retained_target, true_retained_nontarget, rejectedasvshash
	rejectedasvs = tuple(sorted(rejectedasvs))
	return([score] + inputs + retained_vals + estimates + [hash(rejectedasvs)])

def write_specs_and_stats(specs, thresholds, scores, path):
	
	with open(path, "w") as o:
		
		# Write header
		
		scorehead = ["score", "asvs", "target", "nontarget", "rejectedasvs", "retainedasvs", "retained_target", "retained_nontarget", "actual_retainedasvs", "est_true_target", "est_true_nontarget", "est_true_retained_target", "est_true_retained_nontarget, hash_rejectedasvs"]
		
		o.write(",".join([s + "_threshold" for s in specs] + scorehead) + '\n')
		
		# Write lines
		
		for thresh, score in zip(thresholds, scores):
			o.write(",".join(str(v) for v in list(thresh) + score)+'\n')


def get_minimum_thresholds(scores, threshold_combinations, spec):
	
	# Get list of scores only 
	score_list = [ s[0] for s in scores ]
	
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
		
		outtext += "\tPreliminary number of putative exclusions: " + str(scores[min_index][4]) + "\n"
		outtext += "\tProportion of \'verified target\' ASVs retained: " + str(scores[min_index][6]/scores[min_index][2]) + "\n"
		outtext += "\tProportion of \'verified non-target\' ASVs retained: " + str(scores[min_index][7]/scores[min_index][3]) + "\n"
		outtext += "\tFinal number of exclusions: " + str(scores[min_index][8]) + "\n"
		
		n += 1
	
	return(minscore, min_thresholds, outtext)


def output_filtered_haplotypes(counts, min_thresholds, anythreshold, good, bad, file, filename, outdir):
	
	#if(type(min_thresholds[0]) != list):
	#	min_thresholds = [min_thresholds,]
	
	#min_i, thresholds = list(enumerate(min_thresholds))[10]
	
	for min_i, thresholds in enumerate(min_thresholds):
		
		failures = set(itertools.chain.from_iterable([ categorycounting.reject(counts[i], t, anythreshold) for i, t in enumerate(thresholds) ]))
		
		exclude = failures.union(bad) - good
		with open(file) as infasta:
			filen = min_i + 1
			with open(os.path.join(outdir, filename + "_filtered_set" + str(filen) + ".fa"), "w") as outfasta:
				
				for head, seq in SimpleFastaParser(infasta):
					
					if( head not in exclude):
						outfasta.write(">%s\n%s\n" % (head, seq))
