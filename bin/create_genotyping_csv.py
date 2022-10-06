#!/usr/bin/env python3

import os
import re
import csv
import sys

animal = sys.argv[1]

# make list from samples
sam_files = str(sys.argv[2]).split(",")

# create CSV
with open(str(animal) + "_genotyping.csv", 'a') as genotyping_csv:
	
	for i in sam_files:
		# get sample name
		sam_file_basename = os.path.basename(i)
		animal_name = re.sub('.sam', '', sam_file_basename)
	
		with open(i) as tsvfile:
			reader = csv.reader(tsvfile, delimiter='\t')
			for row in reader:
				# get genotype from SAM file
				genotype = row[2]
				
				# write to output CSV
				genotyping_csv.write(animal_name + ',' + genotype + '\n')
