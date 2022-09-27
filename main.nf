#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


workflow {
	
	ch_sample_manifest = Channel
		.fromPath( params.sample_manifest )
		.splitCsv( sep: \t, header: ['bam', 'sample'] )
		.map { row -> tuple( file(row.bam), row.sample )}
	
	CONVERT_BAM_TO_FASTQ ()
	
	ORIENT_FASTQ ()
	
	TRIM_FASTQ ()
	
	INDEX_FASTQ ()
	
	INDEX_GUIDE_FASTA ()
	
	RUN_PBAA ()
	
	RENAME_FOR_LAA ()
	
	CLUSTER_PER_SAMPLE ()
	
	RENAME_CLUSTERS ()
	
	MERGE_PER_ANIMAL_CLUSTERS ()
	
	SHARED_ANIMALS ()
	
	RENAME_PUTATIVE_ALLELE_CLUSTERS ()
	
	PARSE_IPD_GENBANK ()
	
	MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA ()
	
	FILTER_EXACT_GDNA_MATCHES ()
	
	MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE ()
	
	FIND_MUSCLE_CDNA_GDNA_MATCHES ()
	
	RENAME_MUSCLE_CDNA_MATCHES_FASTA ()
	
	EXTRACT_NOVEL_SEQUENCES ()
	
	MERGE_READS ()
	
	CLUSTAL_ALIGN ()
	
	PARSE_DISTANCES ()
	
	CREATE_GENOTYPING_FASTA ()
	
	GENOTYPE_CSS ()
	
	CREATE_GENOTYPING_CSV ()
	
	CREATE_GENOTYPING_PIVOT ()
	
	PRELIMINARY_EXONERATE_PUTATIVE ()
	
	PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE ()
	
	PRELIMINARY_EXONERATE_MERGE_CDS_PUTATIVE ()
	
}

	
process CONVERT_BAM_TO_FASTQ {
	
	// Demultiplexed PacBio CCS reads are provided as unmapped BAM files
	// Convert these CCS files to FASTQ format, renaming FASTQ to use sample names
	// The FASTQ filenames need to contain the locus (dpa, dpb, etc.) in the filenames
	// so the appropriate processing parameters can be applied during orienting, trimming, and pbaa

}

process ORIENT_FASTQ {
	
	// use vsearch orient command to ensure reads are all in the same orientation

}

process TRIM_FASTQ {
	
	// use bbduk to remove 5' and 3' primer sequences
	// failing to remove these sequences leads to artificial alleles with primer sequences at ends
	// need to trim to minlength or else short sequences choke pbaa
	// since sequences are oriented, only need to trim ends in one orientation

}

process INDEX_FASTQ {
	
	// make index from each FASTQ file
	// indexing is necessary for pbaa

}

process INDEX_GUIDE_FASTA {
	
	// make index from each potential guide FASTA file
	// indexing is necessary for pbaa

}

process RUN_PBAA {
	
	// run pbaa clustering on each FASTQ that has been trimmed
	// suppress stderr so you can see which samples fail

}

process RENAME_FOR_LAA {
	
	// pbaa adds extraneous sequence to the output fasta
	// rename to remove

}

process CLUSTER_PER_SAMPLE {
	
	// run bbbmap dedupe.sh on each sample, finding exact matches and absorbing containments.
	// use PacBio clustering parameters from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/ but change edit distance to 0.
	// This should restrict clustering to identical sequences. 
	// outbest parameter is not documented in dedupe.sh but is used in the example from the PacBio clustering and saves one representative sequence per cluster.
	// Since the sequences in the cluster are identical and have an edit distance of 0, this exemplar sequence should be sufficient.
	// 
	// Some samples might not have output. Create empty output file and overwrite if data exists.

}

process RENAME_CLUSTERS {
	
	// prepend sample name to cluster names. This makes it easier to track which sequences are associated with which clusters.

}

process MERGE_PER_ANIMAL_CLUSTERS {
	
	// gymnastics with per-sample merged FASTA to merge them after all have been created - need to wait until rename_clusters completes before running this rule
	// use of lambda function to join all samples into single string in command documented here: https://bitbucket.org/snakemake/snakemake/issues/895/combine-multiple-files-for-input-but
	// this may be broadly useful in cases when we need to merge together files but wait until all the input files exist

}

process SHARED_ANIMALS {
	
	// run bbbmap dedupe.sh on each sample, finding exact matches, absorbing containments, and finding overlaps that are at least 3kb long.
	// this requires a three step process:
	// 1. Run dedupe.sh and output a FASTA file containing all singletons and duplicated sequences
	// 2. Run dedupe.sh again and output a FASTA file containing only singletons.
	// 
	// Now there are singleton sequences in both FASTA files, but sequences that are duplicates (and by definition found in more than one animal)
	// in only one file. So if we find the sequences that are unique *between* these two files, we are left with only the duplicate sequences.
	// So step 3:
	// 3. Remove singleton sequences found in both files
	// 
	// A version of this is described in https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/

}

process RENAME_PUTATIVE_ALLELE_CLUSTERS {
	
	// add integer FASTA id to gdna_match FASTA name
	// this simplifies downstream analyses, where the FASTA description can be complex and the FASTA id is simple
	// since previous steps are nondeterminisitic, add timestamp to id too
	// this makes it easier to clarify which files are used for genotyping

}

process PARSE_IPD_GENBANK {
	
	// use biopython to take IPD Genbank file and separate into FASTA
	// of genomic DNA (intron containing) sequences and cDNA sequences (no intron)

}

process MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA {
	
	// identify putative alleles whose sequences match full-length gDNA sequences already in IPD.
	// save BAM file of reads that map to known IPD alleles and FASTA.gz file of reads that do not.

}

process FILTER_EXACT_GDNA_MATCHES {
	
	// filter mappings to only those that have NM:i:0 (no mismatches)
	// use filterlines.sh tool

}

process MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE {
	
	// among sequences that don't match existing full-length gDNA sequences
	// find ones that extend known cDNA sequences to gDNA versions
	// 
	// rewritten in 23230 to exhaustively pairwise align gDNA sequences to cDNA library with MUSCLE
	// and return alignments where the number of matched nucleotides exactly matches the length of the cDNA
	// this is because the minimap2-based method I used previously failed to find cDNA matches that have a 5nt exon 8
	// 
	// modified in 27309 to minimize search space by first finding cDNA sequences within 50 substitutions of gDNA sequence
	// instead of searching ~1000 cDNA sequences per gDNA sequence, this reduces the search space to ~30-50 in most cases

}

process FIND_MUSCLE_CDNA_GDNA_MATCHES {
	
	// run awk to find gdna sequences that fully match sequence of cdna

}

process RENAME_MUSCLE_CDNA_MATCHES_FASTA {
	
	// input gDNA sequences
	// change the name of the sequence header
	// write output sequences that match cDNAs to new file
	
}

process EXTRACT_NOVEL_SEQUENCES {
	
	// the previous steps processed clusters that perfectly match known IPD cDNA sequences.
	// Because the sequences in this experiment are supported by PacBio cDNA and pbaa gDNA sequences,
	// we can submit them to IPD even if they don't match a current IPD cDNA sequence.
	// 
	// Map gDNA sequences to cDNAs and report reads that don't map

}

process MERGE_READS {
	
	// concatenate files to align with clustal omega

}

process CLUSTAL_ALIGN {
	
	// align reads with clustal omega
	// generate alignment in fasta format and distance matrix
	// distance matrix can be used to parse closest matches to known alleles

}

process PARSE_DISTANCES {
	
	// parse clustal omega distances to find alleles in reference database and cDNA extensions that are nearest neighbors to novel sequences

}

process CREATE_GENOTYPING_FASTA {
	
	// to check accuracy of novel sequences, want to genotype IPD gDNA matches, cDNA matches, and novel sequences
	// need to create FASTA file that contains all of these putative_alleles along with their classification
	// 
	// 27319 - create FASTA file that only has cDNA extensions and novel alleles

}

process GENOTYPE_CSS {
	
	// identify reference sequences that are fully contained in CCS reads
	// use trimmed full-length reads for mapping
	// use higher minratio and smaller maxindel per emails with Brian Bushnell

}

process CREATE_GENOTYPING_CSV {
	
	// create CSV file containing sample name, allele name, and count of reads matching allele

}

process CREATE_GENOTYPING_PIVOT {
	
	// make Excel-formatted pivot table from genotyping CSV

}

process PRELIMINARY_EXONERATE_PUTATIVE {
	
	// run exonerate as in 23188
	// return GFF files with annotations
	// 
	// this preliminary annotation relative to HLA-A is good but often has sequence before the starting methionine and after the stop codon
	// 
	// therefore, parse the annotations with gffread to get protein sequence with correct stop codon. 
	// Use python regexp to trim protein sequence to starting methionine.
	// 
	// then re-run exonerate in protein mode to get more accurate annotations
	// 
	// 27319 - create annotations from FASTA file as named in genotyping file

}

process PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE {
	
	// prepare exonerate GFF output for Geneious
	// need to pass experiment and paths to needed tar.gz files
	// script is in Docker container

}

process PRELIMINARY_EXONERATE_MERGE_CDS_PUTATIVE {
	
	// when Geneious-compatible GFF files are made,
	// the CDS annotations are viewed as independent annotations.
	// This step adds a name and ID annotation to CDS annotations
	// so they view as a single annotation in Geneious and can be automatically translated.
	// 
	// Since we only work with CDS annotations, only print these to final file.

}

