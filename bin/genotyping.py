#!/usr/bin/env python3

import pandas as pd
import sys

animal = sys.argv[1]

# import csv
df = pd.read_csv(str(animal) + "_genotyping.csv"], sep=',', names=['animal', 'genotype'])

# add column to hold counts
df['ct'] = 1

# group by number of times a genotype appears in an animal
df_grouped = df.groupby(['animal', 'genotype'])['ct'].count().reset_index()

# create pivot table
df_pivoted= pd.pivot_table(df_grouped, index=['genotype'], columns=['animal'], values=['ct'])

# export to Excel
# has index column that should be deleted to avoid confusing numbering
df_pivoted.to_excel(snakemake.output[0])
