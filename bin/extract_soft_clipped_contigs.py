import pysam

# read SAM file with pysam and store as python object
samfile = pysam.AlignmentFile(snakemake.input[0], "r")

# write gzip-compressed FASTA
with open(snakemake.output[0], "wb") as f:

	# iterate over reads
	for read in samfile.fetch():
		# only write reads that have a CIGAR mapping
		# and there are soft trims at both ends
		# and aligned length is greater than 
		if read.cigartuples and read.cigartuples[0][0] == 4 and read.cigartuples[-1][0] == 4 and len(read.query_alignment_sequence) >= 2600:
			quality = ''.join(map(lambda x: chr( x+33 ), read.query_alignment_qualities))
			
			# write gzip compressed FASTQ
			# need to set encoding correctly for gzip module to work
			f.write(('@' + read.query_name + '\n').encode("utf-8"))
			f.write((read.query_alignment_sequence + '\n').encode("utf-8"))
			f.write(('+\n').encode("utf-8"))
			f.write((quality + '\n').encode("utf-8"))