import os
import re
import csv

# make list from samples
sam_files = snakemake.params.per_sample_files.split(",")

# create CSV
with open(snakemake.output[0], 'a') as genotyping_csv:
	
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