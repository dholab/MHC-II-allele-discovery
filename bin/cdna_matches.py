#!/usr/bin/env python3

import sys
import pysam
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

mapped_bam = sys.argv[1]
no_gdna_match_fasta = sys.argv[2]

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

# Open a BAM file for reading
bam_file = mapped_bam
bam = pysam.AlignmentFile(bam_file, "rb")

# Read the cDNA reference FASTA
fasta_file = no_gdna_match_fasta

# Open Genbank file
# Save annotated cDNA reference with exon coordinates marked
with open("cdna_matches.gbk", "w") as cdna_extension_handle, open("novel_matches.fasta", "w") as novel_extension_handle:

  # Loop over the reads in the BAM file
  for read in bam:
      # Parse the CIGAR string for the read
      cigar = read.cigar

      # BAM descriptor of 8 is a sequence mismatch
      # retain reads that do not have any mismatches
      mismatch_count = 0
      for op, length in cigar:
          if op == 8:
              mismatch_count = 1

      # if number of mismatches is > 0
      # it is a novel allele
      # output to novel_matches.fasta
      if mismatch_count > 0:
        records = list(SeqIO.parse(fasta_file, "fasta"))

        for record in records:
          # Find matching record
          if read.reference_name == record.name:
            SeqIO.write(record, novel_extension_handle, "fasta")
      
      if mismatch_count == 0:
        
        # Find the corresponding sequence record in the reference FASTA
        records = list(SeqIO.parse(fasta_file, "fasta"))

        for record in records:
          # Find matching record
          if read.reference_name == record.name:

            # Print a summary of exact matches
            print("SPAdes contig %s an exact match to cDNA %s." % (read.reference_name, read.qname))
          
            # Extract the start and stop coordinates of the mapped read
            chrom = read.reference_name
            pos = read.reference_start
            cigar = read.cigartuples

            # Define the exon coordinates as a list of tuples
            exon_coords = []

            # Print the mapped regions for this read
            print("Mapped regions for read %s on %s:" % (read.query_name, chrom))
            gap_start = None
            gap_end = None
            for op, length in cigar:
                if op == 0 or op == 7 or op == 8:
                    # Match or mismatch
                    region_start = pos
                    region_end = pos + length
                    pos = region_end
                    if gap_start is not None:
                        # Print the previous gap
                        print("- %s:%d-%d" % (chrom, gap_start, gap_end))
                        gap_start = None
                        gap_end = None
                    print("- %s:%d-%d" % (chrom, region_start, region_end))
                    
                    # add exon coordinates to list
                    exon_coords.append((region_start, region_end))
                elif op == 1:
                    # Insertion
                    region_start = pos
                    region_end = pos
                    pos += length
                    print("- %s:%d-%d" % (chrom, region_start, region_end))
                elif op == 2:
                    # Deletion
                    if gap_start is None:
                        gap_start = pos
                    gap_end = pos + length
                    pos += length
                elif op == 3:
                    # Skipped region
                    if gap_start is not None:
                        # Print the previous gap
                        print("- %s:%d-%d" % (chrom, gap_start, gap_end))
                        gap_start = None
                        gap_end = None
                    pos += length

            # Create feature annotations for each exon
            exon_features = [SeqFeature(FeatureLocation(start, end), type="exon") for (start, end) in exon_coords]

            # Add the features to the record
            record.features += exon_features

            # Add a source annotation with the matching IPD cDNA sequence
            source_feature = SeqFeature(FeatureLocation(0, len(record.seq)), type="source", qualifiers={"cdna_allele": str(read.query_name)})

            # Add the source feature to the SeqRecord
            record.features.append(source_feature)

            # Set the alphabet to generic DNA
            record.annotations["molecule_type"] = "DNA"

            # Write output file
            SeqIO.write(record, cdna_extension_handle, "genbank")

            # Print the final gap, if there is one
            if gap_start is not None:
                print("- %s:%d-%d" % (chrom, gap_start, gap_end))

  # Close the BAM file
  bam.close()