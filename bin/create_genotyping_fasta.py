import csv
import utils
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import os

# create dictionary of classifications
classified = {}
	
# add cdna matches
with open(snakemake.input[2]) as tsvfile:
	reader = csv.reader(tsvfile, delimiter='\t')
	for row in reader:
		classified[row[0]] = ['extend-cdna', utils.removeSpecialCharacters(row[1])]
		
# add novel matches
if os.stat(snakemake.input[3]).st_size > 0:
	novel_df = pd.read_excel(snakemake.input[3], index_col=0)
	
	for index, row in novel_df.iterrows():
		classified[str(row[0])] = ['novel', utils.removeSpecialCharacters(row[1] + '|' + row[2] + '|' + row[3] + '|' + row[4] + '|' + row[5])]

# create renamed FASTA file with updated names for genotyping
with open(snakemake.output[0], "a") as handle:
	# add IPD gDNA sequences
	with open(snakemake.input[1], "r") as input_handle:
		sequences = SeqIO.parse(input_handle, "fasta")
		SeqIO.write(sequences, handle, "fasta")
	
	# concatenate cDNA extensions and novel sequences with known IPD gDNA sequences
	# this enables genotyping against an expanded gDNA library even when there aren't a huge number of gDNA matches in this specific set of samples		
	for record in SeqIO.parse(snakemake.input[0], "fasta"):
		# get information for matching sequence
		allele_data = classified.get(record.name)
		
		if allele_data is not None:                
			renamed_name = str(allele_data[1] + '_' + allele_data[0] + '_' + record.name)
			renamed_seq = record.seq
			record = SeqRecord(
				  renamed_seq,
				  id=renamed_name,
				  description=''
				  )
				  
			SeqIO.write(record, handle, "fasta")

# create FASTA with only cDNA extension and novel sequences
# Roger prefers having the names of these FASTA sequences matching the genotyping table

with open(snakemake.output[1], "a") as handle:	
	# concatenate cDNA extensions and novel sequences with known IPD gDNA sequences
	# this enables genotyping against an expanded gDNA library even when there aren't a huge number of gDNA matches in this specific set of samples		
	for record in SeqIO.parse(snakemake.input[0], "fasta"):
		# get information for matching sequence
		allele_data = classified.get(record.name)
		
		if allele_data is not None:                
			renamed_name = str(allele_data[1] + '_' + allele_data[0] + '_' + record.name)
			renamed_seq = record.seq
			record = SeqRecord(
				  renamed_seq,
				  id=renamed_name,
				  description=''
				  )
				  
			SeqIO.write(record, handle, "fasta")			
	