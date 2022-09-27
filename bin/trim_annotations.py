import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# import GFF file to dataframe
df = pd.read_csv(snakemake.input[0], sep='\t', names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

# read sequence from FASTA file
with open(snakemake.output[0], "w") as handle:
	for record in SeqIO.parse(snakemake.input[1], "fasta"):

		# filter on single seqid
		seqid_df = df.loc[df['seqid'] == str(record.id)]

		# sort by start coordinate
		seqid_df.sort_values(by=['start'], ascending=True)

		# assuming start of CDS is in-frame
		# determine the frame of each subsequent exon

		# initialize counter to store cumulative nucleotide length
		nuc_length = 0

		# initialize offset to 0 for first annotation
		offset = 0

		# set counter for each CDS annotation for a record
		ct = 0

		for idx, row in seqid_df.iterrows():
			
			# read each GFF row into a variable
			# change these variables when start and stop codon positions change
			seqid = row['seqid']
			source = row['source']
			gff_type = row['type']
			start=row['start']
			end=row['end']
			score=row['score']
			strand=row['strand']
			phase=row['phase']
			attributes=row['attributes']

			# iterate CDS annotation counter
			ct += 1
			
			# write new GFF file
			handle.write(str(seqid) + '\t' 
						+ str(source) + '\t' 
						+ str(gff_type) + '\t' 
						+ str(start) + '\t' 
						+ str(end) + '\t' 
						+ str(score) + '\t'
						+ str(strand) + '\t'
						+ str(phase) + '\t'
						+ str(attributes) + '\n')

			# I tried to handle the stop codon trimming automatically, but could not get accuracy high enough
			# will trim 3' end manually