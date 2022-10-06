#!/usr/bin/env python3

import pandas as pd 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import re
import os
from pathlib import Path
import sys

animal = sys.argv[1]

# test if there is data in novel.fasta
# create output files if empty

if os.stat(str(animal) + "_novel.fasta").st_size == 0:
	Path(str(animal) + "_novel_closest_matches.xlsx").touch()
	Path(str(animal) + "_distances_tmp.txt").touch()
	sys.exit()

# import clustal omega distance file
# remove spaces in sequence identifiers that makes distances.txt file human readable

with open(str(animal) + "_distances.txt", "rt") as fin:
	with open(str(animal) + "_distances_tmp.txt", "wt") as fout:
		for line in fin:
			fout.write(re.sub(' +', ' ', line))

# import into pandas dataframe
df = pd.read_csv(str(animal) + "_distances_tmp.txt", skiprows=1, header=None, sep=" ", index_col=0)

# get row index names
ordered_values = df.index.values.tolist()

# use row_index names as column names
df.columns = ordered_values
df.columns = df.columns.astype(str)

# convert column names that are cdna_matches to cdna names
# this makes it possible to glean lineages for novel alleles that are closest matches to cdna extensions
cdna_dict = { seq_record.name : seq_record.description.split(' ')[1] for seq_record in SeqIO.parse(snakemake.input[2], "fasta")}
df.rename(cdna_dict, inplace = True, axis = 1)

# self-by-self comparisons always have a distance of 0.000000
# for novel alleles, there will never be a perfect match to an existing sequence
# so convert 0.000000	to 1 so self-matches don't report in nsmallest calculation
df.replace(0.000000, 1)

# get names of sequences in novel.fasta
identifiers = [seq_record.id for seq_record in SeqIO.parse(str(animal) + "_reads.fasta", "fasta")]

# create subset df with only novel fasta records as rows
novel_df = df.loc[ identifiers , : ]

# remove columns that are also novel alleles since these are not informative for naming
novel_df.drop(columns = identifiers, inplace=True)

# filter on best three matches
# sort on rows with lowest values (closest matches)
# from https://stackoverflow.com/questions/54923349/top-3-values-per-row-in-pandas
c = ['Self', 'Closest','2','3', '4', '5']
ranked_df = (novel_df.apply(lambda x: pd.Series(x.nsmallest(6).index, index=c), axis=1)
		.reset_index())

# drop 'self' match which will always be the best match
top5_df = ranked_df.drop(columns = ['Self'])
top5_df.columns = top5_df.columns.astype(str)

# export to excel
# top5_df.to_excel('novel.xlsx')

# iterate over each entity in top5_df
# add distance in parentheses after each
# the logic took me a while to figure out
# first lookup each coordinate in top5_df
# then find the corresponding value in novel_df
# and append it as a parenthetical

for index, row in top5_df.iterrows():
	# update value of each combination
	top5_df.at[index, 'Closest'] = top5_df.at[index, 'Closest'] + ' (' + str(novel_df.at[row['0'], row['Closest']].round(3)) + ')'
	top5_df.at[index, '2'] = top5_df.at[index, '2'] + ' (' + str(novel_df.at[row['0'], row['2']].round(3)) + ')'
	top5_df.at[index, '3'] = top5_df.at[index, '3'] + ' (' + str(novel_df.at[row['0'], row['3']].round(3)) + ')'
	top5_df.at[index, '4'] = top5_df.at[index, '4'] + ' (' + str(novel_df.at[row['0'], row['4']].round(3)) + ')'
	top5_df.at[index, '5'] = top5_df.at[index, '5'] + ' (' + str(novel_df.at[row['0'], row['5']].round(3)) + ')'

top5_df.to_excel(str(animal) + "_novel_closest_matches.xlsx")
