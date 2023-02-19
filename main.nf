#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// input channels
	ch_sample_manifest = Channel
		.fromPath( params.sample_manifest )
		.splitCsv( header: ['fastq', 'sample'] )
		.map { row -> tuple( "${params.spades_contig_folder}/${row.fastq}", row.sample) }
	
	// Workflow steps
	
	PREFIX_FASTA (
		ch_sample_manifest
	)

	ORIENT_FASTA (
		PREFIX_FASTA.out.prefixed_fasta
	)

	MAP_FASTA (
		ORIENT_FASTA.out.oriented_fasta
	)

	TRIM_TO_PACBIO_AMPLICONS (
		MAP_FASTA.out.mapped_bam
	)

	FILTER_HARD_CLIPPED_AMPLICONS (
		TRIM_TO_PACBIO_AMPLICONS.out.trimmed_bam
	)
	
	// RENAME_CLUSTERS (
	// 	FILTER_HARD_CLIPPED_AMPLICONS.out.clipped_fasta
	// )
	
	MERGE_PER_SAMPLE_CLUSTERS (
		FILTER_HARD_CLIPPED_AMPLICONS.out.clipped_fasta
			.map { fasta, sample-> fasta }
			.collect()
	)
	
	SHARED_ANIMALS (
		MERGE_PER_SAMPLE_CLUSTERS.out.merged_fasta
	)
	
	RENAME_PUTATIVE_ALLELE_CLUSTERS (
		SHARED_ANIMALS.out.putative
	)
	
	CLASSIFY_PUTATIVE (
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out.renamed_clusters,
		params.classify_genbank
	)

	GENOTYPE_PUTATIVE(
		PREFIX_FASTA.out.prefixed_fasta
			.map { fasta, sample-> fasta }
			.collect(),
		CLASSIFY_PUTATIVE.out.all_fasta,
		
	)
	
	// MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA (
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out.renamed_clusters,
	// 	PARSE_IPD_GENBANK.out.ipd_gdna
	// )
	
	// FILTER_EXACT_GDNA_MATCHES (
	// 	MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out.all_mappings,
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out.renamed_clusters,
	// )
	
	// DEFINE_CDNA_MATCHES_AND_NOVELS (
	// 	FILTER_EXACT_GDNA_MATCHES.out.no_gdna_matches_fasta,
	// 	PARSE_IPD_GENBANK.out.ipd_cdna
	// )
	
	// MERGE_CDNA_MATCHES_AND_NOVELS (
	// 	PARSE_IPD_GENBANK.out.ipd_gdna,
	// 	DEFINE_CDNA_MATCHES_AND_NOVELS.out.cdna_extensions_fasta,
	// 	DEFINE_CDNA_MATCHES_AND_NOVELS.out.novel_alleles_fasta
	// )
	
	// CLUSTAL_ALIGN (
	// 	MERGE_CDNA_MATCHES_AND_NOVELS.out.merged_for_clustal_fasta
	// )
	
	// FIND_CLOSEST_MATCHES (
	// 	CLUSTAL_ALIGN.out.clustal_distances,
	// 	DEFINE_CDNA_MATCHES_AND_NOVELS.out.novel_alleles_fasta,
	// 	DEFINE_CDNA_MATCHES_AND_NOVELS.out.cdna_extensions_fasta
	// )
	
	// CREATE_GENOTYPING_FASTA (
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out.renamed_clusters,
	// 	PARSE_IPD_GENBANK.out.ipd_gdna,
	// 	DEFINE_CDNA_MATCHES_AND_NOVELS.out.cdna_extensions_fasta,
	// 	FIND_CLOSEST_MATCHES.out.closest_matches
	// )
	
	// PRELIMINARY_EXONERATE_PUTATIVE (
	// 	CREATE_GENOTYPING_FASTA.out.new_allele
	// )
	
	// PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE (
	// 	PRELIMINARY_EXONERATE_PUTATIVE.out.gff
	// )
	
	// PRELIMINARY_EXONERATE_MERGE_CDS_PUTATIVE (
	// 	PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE.out.processed_gff
	// )
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Derivative parameters, mostly for making specific results folders
// params.orient_fastq = params.results + "/" + "00-orient-fastq"
// params.map_fastq = params.results + "/" + "01-map-fastq"
// params.trim_to_pacbio_amplicons = params.results + "/" + "02-trim-to-pacbio-amplicons"
// params.filter_hard_clipped_amplicons = params.results + "/" + "03-filter-hard-clipped-amplicons"
// params.sample_clusters = params.results + "/" + "04-cluster_per_sample"
// params.merged_clusters = params.results + "/" + "05-merged_clusters"
// params.shared_clusters = params.results + "/" + "06-shared_clusters"
// params.ipd_ref_sep = params.results + "/" + "ipd_ref_separate"
// params.gdna_identical = params.results + "/" + "07-gdna-identical"
// params.cdna_identical = params.results + "/" + "08-muscle-cdna-identical"
// params.novel_alleles = params.results + "/" + "09-novel"
// params.genotyping = params.results + "/" + "10-genotyping"
// params.putative_new = params.results + "/" + "11-putative-new-alleles"
// params.ipd_refs = params.classify_resources + "/" + "*.gbk"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process PREFIX_FASTA {
	/*
	* Process: PREFIX_FASTA
	* Description: Prefix FASTA reads with sample name
	* This makes it easier to track reads throughout the workflow
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	
	// Set process settings for resource allocation and error handling
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	// Define input parameters for the process
	input:
	tuple path(fasta), val(sample)
	
	// Define output parameters for the process
	output:
	tuple path("${sample}.renamed.fasta"), val(sample), emit: prefixed_fasta
	
	// Define the script to execute for the process
	script:
	"""
	rename.sh in=${fasta} prefix=${sample} addprefix=t out="${sample}.renamed.fasta"
	"""
}

process ORIENT_FASTA {
	/*
	* Process: ORIENT_FASTA
	* Description: Use vsearch orient command to ensure reads are all in the same orientation.
	* This avoids complex reverse complementing some contigs in subsequent steps.
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	
	// Set process settings for resource allocation and error handling
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	// Define input parameters for the process
	input:
	tuple path(fasta), val(sample)
	
	// Define output parameters for the process
	output:
	tuple path("${sample}.fasta"), val(sample), emit: oriented_fasta
	
	// Define the script to execute for the process
	script:
	"""
	vsearch --orient ${fasta} \
	--db ${params.combined_reference} \
	--fastaout ${sample}.fasta
	"""
}

process MAP_FASTA {
	
	/*
	* Process: MAP_FASTA
	* Description: Use minimap2 to map SPAdes contigs to combined DPA, DPB, DQA, DQB, DRB reference FASTA.
	* This reference sequence contains 1kb of flanking sequence at end of each gene,
	* which is needed to ensure that indels near the PacBio amplicon ends map correctly.
	* After mapping, use samtools to sort the minimap2 output and convert to BAM file.
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	
	// Set process settings for resource allocation and error handling
	cpus 1
	memory '2.5 GB'
	errorStrategy 'retry'
	maxRetries 4
	
	// Define input parameters for the process
	input:
	tuple path(fasta), val(sample)
	
	// Define output parameters for the process
	output:
	tuple path("${sample}_mapped.bam"), val(sample), emit: mapped_bam
	
	// Define the script to execute for the process
	script:
	"""
	minimap2 -ax asm20 \
	${params.combined_reference} \
	${fasta} \
	--sam-hit-only --eqx \
	| samtools sort - \
	> ${sample}_mapped.bam
	"""
}

process TRIM_TO_PACBIO_AMPLICONS {
	/*
	* Process: TRIM_TO_PACBIO_AMPLICONS
	* Description: Use samtools ampliconclip to trim mapped sequences to PacBio amplicon coordinates.
	* This requires a BED file with the amplicon coordinates relative to the 1kb extended reference sequences.
	* When using ampliconclip, it is necessary to take advantage of strand information to trim the correct end of the sequence.
	* A 9100bp tolerance is required because this is where the PacBio amplicon is relative to the 1kb extended DRB reference sequence.
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	
	// Set process settings for resource allocation and error handling
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	// Define input parameters for the process
	input:
	tuple path(bam), val(sample)
	
	// Define output parameters for the process
	output:
	tuple path("${sample}_trimmed.bam"), val(sample), emit: trimmed_bam
	
	// Define the script to execute for the process
	script:
	"""
	samtools ampliconclip \
	-b ${params.primer_bed} \
	--hard-clip \
	--both-ends \
	${sample}_mapped.bam \
	--tolerance 9100 \
	--strand \
	--clipped \
	| samtools sort \
	> ${sample}_trimmed.bam \
	&& samtools index ${sample}_trimmed.bam
	"""
	
}

process FILTER_HARD_CLIPPED_AMPLICONS {
	/*
	* FILTER_HARD_CLIPPED_AMPLICONS
	*
	* This process filters sequences that span the entire PacBio amplicon.
	* After running ampliconclip, sequences that have hard clips at both ends are the ones we want.
	* This process uses a custom script to filter these sequences.
	*/

	
	/*
	 * Other parameters:
	 *   - tag: tag for the process
	 */
	tag "${sample}"

	// Set process settings for resource allocation and error handling
	cpus 1
	errorStrategy 'retry'
	maxRetries 4

	/*
	 * Input parameters:
	 *   - bam: path to the BAM file containing mapped sequences hard clipped to the PacBio amplicon
	 *   - sample: sample name
	 */
	input:
	tuple path(bam), val(sample)
	
	/*
	 * Output parameters:
	 *   - fasta: path to the output FASTA file containing filtered sequences with hard clipped sequences removed
	 *   - sample: sample name
	 */
	output:
	tuple path("${sample}_filtered.fasta"), val(sample), emit: clipped_fasta
	

	/*
	 * Execution script:
	 *   - Indexes the trimmed BAM file so it can be parsed by pysam
	 *   - Calls a Python script to filter sequences with hard clips at both ends
	 *   - Converts the filtered BAM file to a FASTA file
	 */

	script:
	"""
	samtools index ${sample}_trimmed.bam
	
	python ${baseDir}/bin/filter_hard_clipped_ends.py ${sample}_trimmed.bam ${sample}_filtered.bam

	reformat.sh in=${sample}_filtered.bam out=${sample}_filtered.fasta
	"""
	
}

process RENAME_CLUSTERS {
	// Prepend sample name to cluster names. This makes it easier to track which sequences are associated with which clusters.

	tag "${sample}"
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(sample)
	
	output:
	tuple path("${sample}.fasta.gz"), val(sample), emit: renamed_fasta
	
	script:
	"""
	rename.sh -Xmx1g \
	in=${fasta} \
	out=${sample}.fasta.gz \
	prefix=${sample} \
	addprefix=t \
	threads=${task.cpus}
	"""
}

process MERGE_PER_SAMPLE_CLUSTERS {
	/*
	This process merges per-sample clusters into a single file.
	The input to this process is a list of per-sample merged FASTA files.
	It uses zcat and gzip to merge all the files and outputs a single FASTA file containing all the merged clusters.
	The output file is published to the `merged_clusters` directory specified in the parameters.
	*/
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(mamu_list)
	
	output:
	path("merged_clusters.fasta.gz"), emit: merged_fasta
	
	script:
	"""
	cat *.fasta | gzip > merged_clusters.fasta.gz
	"""
}

process SHARED_ANIMALS {
	/*
	* Run bbbmap dedupe.sh on each sample, finding exact matches, absorbing containments, and finding overlaps that are at least 3kb long.
	* This requires a three-step process:
	* 1. Run dedupe.sh and output a FASTA file containing all singletons and duplicated sequences
	* 2. Run dedupe.sh again and output a FASTA file containing only singletons.
	* 
	* Now there are singleton sequences in both FASTA files, but sequences that are duplicates (and by definition found in more than one animal)
	* in only one file. So if we find the sequences that are unique *between* these two files, we are left with only the duplicate sequences.
	* So step 3:
	* 3. Remove singleton sequences found in both files
	* 
	* A version of this is described in https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/
	*/
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(fasta)
	
	output:
	path("putative_alleles_temp.fasta"), emit: putative
	
	script:
	"""
	dedupe.sh -Xmx1g \
	in=${fasta} \
	outbest=all.fasta \
	am=t ac=f arc=t fo c fcc nam=4 threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=${fasta} \
	out=unique.fasta \
	am=t ac=f arc=t fo fcc uniqueonly=t threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=all.fasta,unique.fasta \
	out=shared.fasta \
	ac=f uniqueonly=t threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=shared.fasta \
	out=putative_alleles_temp.fasta \
	ac=t threads=${task.cpus}
	"""
	
}

process RENAME_PUTATIVE_ALLELE_CLUSTERS {
	
	// Add integer FASTA id to gdna_match FASTA name to simplify downstream analyses.
	// A timestamp is also added to the id to help identify the files used for genotyping.
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(putative)
	
	output:
	path("putative_alleles.fasta"), emit: renamed_clusters
	
	script:
	"""
	#!/usr/bin/env python3
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	from datetime import datetime
	
	with open("putative_alleles.fasta", "w") as handle:
		for idx, record in enumerate(SeqIO.parse("${putative}", "fasta")):
			now = datetime.now()
			time = now.strftime("%Y%m%d%H%M%S")
			record.id = str(time) + '-' + str(idx)
			SeqIO.write(record, handle, "fasta")
	"""

}

process CLASSIFY_PUTATIVE {
	
	// Map known IPD sequences to putative alleles and classify as gDNA (exact or extension), 
	// cDNA, novel, or mystery (no match to any known IPD sequence)
	// Output Genbank files for each classification as well as a merged file with all sequences
	// Also output a FASTA file containing all of the putative alleles with informative names
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(putative_alleles)
	path(ipd_genbank)
	
	output:
	path("gdna_ipd.gbk"), emit: gdna_ipd
	path("gdna_extend.gbk"), emit: gdna_extend
	path("cdna_extend.gbk"), emit: cdna_extend
	path("novel.gbk"), emit: novel
	path("mystery.gbk"), emit: mystery
	path("all.gbk"), emit: all_genbank
	path("all.fasta"), emit: all_fasta
	path("ipd.fasta"), emit: ipd_fasta

	script:
	"""
	# convert IPD Genbank file to FASTA
	ipd_genbank_to_fasta.py ${ipd_genbank}

	# map putative alleles to IPD sequences
	minimap2 \
	-ax splice ${putative_alleles} \
	ipd.fasta \
	--eqx \
	--sam-hit-only \
	| samtools sort > mapped.bam \
	&& samtools index mapped.bam

	# create empty output files
	# these will be appended to in the next step
	touch gdna_ipd.gbk gdna_extend.gbk cdna_extend.gbk novel.gbk mystery.gbk all.gbk all.fasta

	# classify putative alleles
	classify_putative.py ${putative_alleles}
	"""

}

process GENOTYPE_PUTATIVE {
	
	// Concatenate IPD database and putative allele FASTA files with nearest IPD match
	// Map each animal's clusters to the concatenated FASTA file
	// Filter mappings to only include reads that have an exact match to gDNA sequences (no mismatches).
	// Use the `filterlines.sh` tool to filter SAM files by NM:i:0 tag and `filterbyname.sh` to filter out
	// FASTA sequences that have matches in the filtered SAM files.
	
	cpus 40
	errorStrategy 'retry'
	maxRetries 4
	
	input:
		path(original_renamed_fasta)
		path(all_putative_fasta)
		
	output:
		path("genotypes.sam"), emit: genotypes_sam
		path("genotypes.csv"), emit: csv
	
	script:
	"""
	# replace sequence names in all_fasta with sequence descriptions
	create_renamed_fasta.py ${all_putative_fasta}
	
	# merge renamed source FASTA into single file
	# retain only animal name in sequence headers
	cat *.renamed.fasta | sed 's/_.*//' > spades.fasta

	# map all clusters to all.fasta
	minimap2 -t ${task.cpus} \
	-ax asm20 all_putative_fasta_renamed.fasta \
	spades.fasta \
	--eqx \
	--sam-hit-only \
	| samtools sort > mapped.bam \
	&& samtools index mapped.bam

	# filter SAM file to only include reads that have an exact match to gDNA sequences (no mismatches)
	filterlines.sh \
		in=mapped.bam \
		out=genotypes.sam \
		names=NM:i:0 \
		substring=t \
		include=t

	# create output genotyping CSV file that can be converted to an Excel Pivot Table
	echo "sample\tallele" > genotypes.csv
	cut -f 1,3 genotypes.sam \
	| awk '{OFS="," ; print \$1, \$2}' \
	>> genotypes.csv
	"""
	
}

process DEFINE_CDNA_MATCHES_AND_NOVELS {

	// Run Python script to define cDNA matches and novel alleles
	// Input files are FASTA files of IPD cDNA sequences and putative allele sequences that did not match existing IPD genomic DNA sequence
	// Output files are two FASTA files: one containing cDNA matches and one containing novel cDNA alleles

	cpus 1
	errorStrategy 'retry'
	maxRetries 4

	input:
	path(no_gdna_matches)
	path(ipd_cdna)
	
	output:
	path("cdna.fasta"), emit: cdna_extensions_fasta
	path("novel.fasta"), emit: novel_alleles_fasta
	
	script:
	"""
	# Use minimap2 and samtools to generate a BAM file of mapped cDNA reads
	minimap2 -ax splice ${no_gdna_matches} ${ipd_cdna} --eqx --sam-hit-only | samtools sort > mapped.bam && samtools index mapped.bam

	# Run a Python script to determine which cDNA sequences are matches and which are novel alleles
	cdna_matches.py mapped.bam ${no_gdna_matches}
	"""	
}

process MERGE_CDNA_MATCHES_AND_NOVELS {

	/*
	Concatenate the input files to prepare for alignment with Clustal Omega.
	*/

	errorStrategy 'retry'
	maxRetries 4

	input:
	path(gdna_ref)
	path(cdna_match)
	path(novel)

	output:
	path("merged_reads.fasta"), emit: merged_for_clustal_fasta
	
	script:
	"""
	cat ${gdna_ref} ${cdna_match} ${novel} > merged_reads.fasta
	"""

}

process CLUSTAL_ALIGN {

	// This process aligns reads with Clustal Omega, generates an alignment in FASTA format, and produces a distance matrix.
	// The distance matrix can be used to identify the closest matches to known alleles.
	
	cpus 40
	executor 'local'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(merged)
	
	output:
	path("aligned.fasta"), emit: clustal_aligned
	path("distances.txt"), emit: clustal_distances
	
	script:
	"""
	clustalo \
	--infile=${merged} \
	--outfile=aligned.fasta \
	--distmat-out=distances.txt \
	--threads=${task.cpus} \
	--full
	"""
}

process FIND_CLOSEST_MATCHES {
	
	// Find the closest matches in the reference database and cDNA extensions to novel sequences using Clustal Omega distances
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(distances)
	path(novel)
	path(cdna_matches)
	
	output:
	path("novel_closest_matches.xlsx"), emit: closest_matches
	path("distances_tmp.txt"), emit: distances
	
	script:
	"""
	parse_clustalo_distances.py ${novel} ${cdna_matches}
	"""

}

process CREATE_GENOTYPING_FASTA {

    /*
    Create a FASTA file that contains putative alleles along with their classification 
    for genotyping IPD gDNA matches, cDNA matches, and novel sequences. 
    The file will contain only cDNA extensions and novel alleles.
    */
	
    errorStrategy 'retry'
    maxRetries 4
	
    input:
    path(putative)
    path(gdna_ref)
    path(cdna_matches)
    path(closest_matches)
	
    output:
    path("classified.fasta"), emit: cdna_extension_fasta
    path("putative.fasta"), emit: putative_novel_allele_fasta
	
    script:
    """
    create_genotyping_fasta.py ${putative} ${gdna_ref} ${cdna_matches}
    """
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

	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(animal)
	
	output:
	tuple path("*mapped.gff"), val(animal), emit: gff
	tuple path("*gdna_single_temp.fasta"), val(animal)
	
	script:
	"""
	#!/usr/bin/env python3
	
	import os
	import subprocess
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	import csv
	
	# create output file in case it is empty
	if os.stat("${fasta}").st_size == 0:
		subprocess.run('touch ${fasta}', shell=True)
	else:
		# read input FASTA line-by-line
		for record in SeqIO.parse("${fasta}", "fasta"):
			if 'dpa' in (record.name).lower():
				mrna_reference = "${params.dpa_mrna_reference}"
				cds_annotation = "${params.dpa_cds_annotation}"
			elif 'dpb' in (record.name).lower():
				mrna_reference = "${params.dpb_mrna_reference}"
				cds_annotation = "${params.dpb_cds_annotation}"
			elif 'dqa' in (record.name).lower():
				mrna_reference = "${params.dqa_mrna_reference}"
				cds_annotation = "${params.dqa_cds_annotation}"
			elif 'dqb' in (record.name).lower():
				mrna_reference = "${params.dqb_mrna_reference}"
				cds_annotation = "${params.dqb_cds_annotation}"
			elif 'drb' in (record.name).lower():
				mrna_reference = "${params.drb_mrna_reference}"
				cds_annotation = "${params.drb_cds_annotation}"
			
			with open("${animal}" + "_gdna_single_temp.fasta", "w") as output_handle:
				SeqIO.write(record, output_handle, "fasta")
			
			# run exonerate        
			subprocess.run('exonerate \
			--showtargetgff \
			--showalignment FALSE \
			--showvulgar FALSE \
			--model cdna2genome \
			--query ' + mrna_reference + ' \
			--target ${animal}_gdna_single_temp.fasta \
			--refine full \
			--annotation ' + cds_annotation + ' \
			>> ${animal}_mapped.gff', shell=True)
	
	"""

}

process PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE {
	
	// prepare exonerate GFF output for Geneious
	// need to pass experiment and paths to needed tar.gz files
	// script is in Docker container
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(gff), val(animal)
	
	output:
	tuple path("*processed.gff"), val(animal), emit: processed_gff
	
	script:
	"""
	bash /scripts/process-gff.sh \
	-e ${gff}  \
	-p 21295-exonerate_gff_to_alignment_gff3.pl \
	-o ${animal}_processed.gff
	"""

}

process PRELIMINARY_EXONERATE_MERGE_CDS_PUTATIVE {
	
	// when Geneious-compatible GFF files are made,
	// the CDS annotations are viewed as independent annotations.
	// This step adds a name and ID annotation to CDS annotations
	// so they view as a single annotation in Geneious and can be automatically translated.
	// 
	// Since we only work with CDS annotations, only print these to final file.
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(gff), val(animal)
	
	output:
	path("*putative.gff")
	
	shell:
	'''
	awk \'{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' !{gff} >> !{animal}_putative.gff
	'''

}
// --------------------------------------------------------------- //


