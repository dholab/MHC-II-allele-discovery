#!/usr/bin/env python3

import sys
import pysam
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

# Open a BAM file for reading
bam_file = sys.argv[1]
bam = pysam.AlignmentFile(bam_file, "rb")

# Make set of sequences that match known cDNA sequences
cdna_matches = set()

# Read the cDNA reference FASTA
fasta_file = sys.argv[2]

# Open Genbank files for cDNA extensions and novel sequences
with open("cdna.fasta", "w") as cdna_output_handle, open("novel.fasta", "w") as novel_output_handle:

  # Loop over the reads in the BAM file
  for read in bam:
      
      # Skip unmapped reads
      if read.is_unmapped:
        continue

      # Parse the CIGAR string for the read
      cigar = read.cigar

      # BAM descriptor of 8 is a sequence mismatch
      # retain reads that do not have any mismatches
      mismatch_count = 0
      for op, length in cigar:
          if op == 8:
              mismatch_count = 1

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
            SeqIO.write(record, cdna_output_handle, "fasta")

            # Add sequence name to matched cDNA set
            cdna_matches.add(read.reference_name)

  # determine which reference sequences do not have cdna_matches 
  bam_ref_names = set(bam.references)
  novel_references = bam_ref_names.difference(cdna_matches)
  
  # Use SeqIO.parse() to read the FASTA file
  records = list(SeqIO.parse(fasta_file, "fasta"))
  for record in records:
      # Check if the record name matches the desired name
      if record.id in novel_references:
          # Print novel output file
          SeqIO.write(record, novel_output_handle, "fasta")

  # Close the BAM file
  bam.close()