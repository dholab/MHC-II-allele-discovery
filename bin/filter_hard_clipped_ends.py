#! /usr/bin/env python3

import pysam
import argparse

# create an ArgumentParser object
parser = argparse.ArgumentParser(description="Parser for filter_hard_clipped_ends.py")

# add an argument for the input file
parser.add_argument("input_file", help="samtools ampliconclip output BAM file")

# add an argument for the output file
parser.add_argument("output_file", help="output BAM file")

# parse the command-line arguments
args = parser.parse_args()

input_bam = pysam.AlignmentFile(args.input_file, "rb")  # open BAM file for reading

output_bam = pysam.AlignmentFile(args.output_file, "wb", template=input_bam)

for read in input_bam:
    
    # intialize counter
    hard_clip_count = 0
    
    if read.cigartuples is not None:
      # iterate over cigartuples
      for cigar_op, cigar_len in read.cigartuples:

        # cigar operation for hard clipping is 5
        if cigar_op == 5:
          hard_clip_count += 1

    # print contigs where there are hard clips at both ends
    if hard_clip_count == 2:
        output_bam.write(read)

input_bam.close()  # close BAM file
output_bam.close()