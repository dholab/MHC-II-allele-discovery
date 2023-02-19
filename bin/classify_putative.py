#!/usr/bin/env python3

import pysam
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
import sys

# command line parameter of putative alleles
sequences_to_classify = sys.argv[1]

def removeSpecialCharacters(in_str, special_characters='*|: ', replace_character='_'):
    '''remove specified special characters from input str'''
    
    out_str = in_str.translate ({ord(c): replace_character for c in special_characters})
    
    return out_str

def annotate_putative(ipd_mapping, putative_allele, match_type, mismatch_count):
    '''
    annotate putative allele with information from best ipd_mapping
    '''
    # Extract the start and stop coordinates of the mapped read
    chrom = ipd_mapping.reference_name
    pos = ipd_mapping.reference_start
    cigar = ipd_mapping.cigartuples

    # Define the exon coordinates as a list of tuples
    exon_coords = []

    gap_start = None
    gap_end = None
    for op, length in cigar:
        if op == 0 or op == 7 or op == 8:
            # Match or mismatch
            region_start = pos
            region_end = pos + length
            pos = region_end
            if gap_start is not None:
                gap_start = None
                gap_end = None
            
            # add exon coordinates to list
            exon_coords.append((region_start, region_end))
        elif op == 1:
            # Insertion
            region_start = pos
            region_end = pos
            pos += length
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
                gap_start = None
                gap_end = None
            pos += length

    # Create a Seq object from the sequence string
    annotated_allele = SeqRecord(putative_allele.seq, id=putative_allele.name, description=putative_allele.description)

    # Set the alphabet to generic DNA
    annotated_allele.annotations["molecule_type"] = "DNA"

    # Create feature annotations for each exon
    exon_features = [SeqFeature(FeatureLocation(start, end), type="match") for (start, end) in exon_coords]

    # Add the features to the record
    annotated_allele.features += exon_features

    # Add a source annotation with the closest IPD sequence
    source_feature = SeqFeature(FeatureLocation(0, len(putative_allele.seq)), type="source", qualifiers={"Closest IPD allele": str(ipd_mapping.query_name), "Mismatches from closest IPD allele": str(mismatch_count), "Category of putative allele": match_type})

    # Add the source feature to the SeqRecord
    annotated_allele.features.append(source_feature)

    # Set the record name to the reference name
    annotated_allele.id = str(ipd_mapping.reference_name)
    annotated_allele.name = str(ipd_mapping.reference_name)
    annotated_allele.description = str(ipd_mapping.query_name) + '|' + str(ipd_mapping.reference_name)

    return(annotated_allele)

# Read the putative allele FASTA from the command line
fasta_file = sequences_to_classify

# Find the corresponding sequence record in the reference FASTA
records = list(SeqIO.parse(fasta_file, "fasta"))

# Open a BAM file for reading
bam_file = "mapped.bam"
bam = pysam.AlignmentFile(bam_file, "rb")

# Interrogate each reference sequence (= putative allele) individually
# This allows the best scoring mapping for our purposes to be retained

reference_names = bam.references
print(reference_names)

# Iterate over each reference name in the BAM file
# Put each one into its own loop
# That way we can break out of the loop when we have a gDNA match
# Before all the mappings have been examined, since the mapping can't be better
# Than an allele we've already described

for record in records:
    print(record.name)

    ## STEP 1: Filter BAM file on best matches for each putative allele ##

    # initialize variables used to rank outputs
    is_gdna = False
    is_gdna_extension = False
    is_cdna_extension = False
    is_novel = False

    # since pysam parses BAM as an iterator
    # need to create new bam object for each iteration
    bam = pysam.AlignmentFile(bam_file, "rb")

    # track best mapping for each putative allele
    # set to an unreasonably high default value
    # so all mismatches from reads will be lower
    best_match = 999999

    for read in bam:

        # find mappings for each putative_allele
        if read.reference_name == record.name:
            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Parse the CIGAR string for the read
            cigar = read.cigar

            # BAM descriptor of 8 is a sequence mismatch
            # retain reads that do not have any mismatches
            # these could be either matches to existing gDNA sequences already in IPD
            # these will not have any deletions relative to the reference
            # or they will be cDNA extensions that have deletions relative to the reference
            # at intron locations

            # also test if gDNA match or cDNA match
            # skipped region in reference is cigar description 3
            mismatch_count = 0
            skipped_count = 0
            soft_clipped_count = 0

            for op, length in cigar:
                if op == 8:
                    mismatch_count += 1

            for op, length in cigar:
                if op == 3:
                    skipped_count += 1

            for op, length in cigar:
                if op == 4:
                    soft_clipped_count += 1

            # use counts from above to classify putative alleles
            # first case is if the putative_allele is fully contained in an IPD gDNA sequence
            # with no soft clipping at ends of 
            # this putative allele cannot have a better match
            if (len(record.seq) <= read.infer_read_length() 
                and mismatch_count == 0 
                and skipped_count == 0):
                is_gdna = True
                best_match = mismatch_count
                best_read = read

            # if the putative allele matches an existing IPD gDNA sequence perfectly
            # but the putative allele sequence is longer than the IPD gDNA sequence
            # the putative allele sequence extends the existing gdna sequence
            elif (len(record.seq) > read.infer_read_length() 
                and mismatch_count == 0 
                and skipped_count == 0):
                is_gdna_extension = True
                best_match = 0
                best_read = read

            # if the putative allele has no mismatches but the reference has skipped regions
            # the putative allele extends an existing IPD cDNA sequence to a gDNA sequence
            elif (len(record.seq) > read.infer_read_length() 
            and mismatch_count == 0):
                is_cdna_extension = True
                best_match = mismatch_count
                best_read = read

            # next find putative alleles that have mismatches relative to IPD
            # these are putative novel alleles
            # we want to save the read that has the lowest number of mismatches relative to an existing IPD sequence
            else:
                if mismatch_count < best_match:
                    best_match = mismatch_count
                    best_read = read
                    is_novel = True
    
    ## STEP 2: Parse best matches for each putative allele ##

    # if the putative allele sequence is already in IPD and cannot be improved
    if is_gdna:
        annotated_seq = annotate_putative(best_read, record, 'gdna_ipd', best_match)
        # write outputs
        with open("gdna_ipd.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.fasta", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "fasta")
        continue

    # if the putative allele sequence is already in IPD and can be improved
    if is_gdna_extension:
        annotated_seq = annotate_putative(best_read, record, 'gdna_extend', best_match)
        with open("gdna_extend.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.fasta", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "fasta")
        continue

    # if the putative allele sequence is in IPD as a cDNA sequence and can be extended
    if is_cdna_extension:
        annotated_seq = annotate_putative(best_read, record, 'cdna_extend', best_match)
        # write outputs
        with open("cdna_extend.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.fasta", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "fasta")
        continue
    
    # if the putative allele sequence is a novel allele
    # report the closest IPD match
    if is_novel:
        annotated_seq = annotate_putative(best_read, record, 'novel', best_match)
        # write outputs
        with open("novel.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.gbk", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "genbank")
        with open("all.fasta", "a") as output_handle:
            SeqIO.write(annotated_seq, output_handle, "fasta")
        continue

    # handle mystery putative alleles that are not classified
    if is_gdna == False and is_gdna_extension == False and is_cdna_extension == False and is_novel == False:
        # mystery sequence that doesn't match anything in IPD
        # Print mystery output file
        with open("mystery.fasta", "a") as output_handle:
            SeqIO.write(record, output_handle, "fasta")