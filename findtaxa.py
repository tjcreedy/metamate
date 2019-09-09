#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Function for reading and checking taxon file"""

# Imports

from io import StringIO
import os
import csv
import sys
from collections import defaultdict


# Global variables

# Class definitions

# Function definitions

def parse_taxa(taxafile, names):
	
	taxa = defaultdict(set)
	allnames = set()
	
	# Load in taxon file csv to dict
	with open(taxafile, 'r') as fh:
		reader = csv.reader(fh)
		for row in reader:
			taxa[row[1]].add(row[0])
			allnames.add(row[0])
	
	# Check that all of the names are represented
	if(set(names) != allnames):
		sys.exit("Error: haplotype names in taxon file do not completely match haplotype names in zotu file")
	
	return(taxa)

