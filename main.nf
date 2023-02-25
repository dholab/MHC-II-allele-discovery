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
		FILTER_HARD_CLIPPED_AMPLICONS.out.clipped_fasta
			.map { fasta, sample-> fasta }
			.collect(),
		CLASSIFY_PUTATIVE.out.all_fasta,
		
	)
}
// --------------------------------------------------------------- //


// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process PREFIX_FASTA {
	/*
	* Process: PREFIX_FASTA
	* Description: Prefix SPAdes contigs in FASTA format with sample name followed by an underscore
	* e.g., NODE_2611_length_3380_cov_80.013210 -> r12083_NODE_2611_length_3380_cov_80.013210
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
	* Description: Use vsearch orient command to ensure reads are all in the same orientation relative to MHC class II reference exemplars.
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
	* Save mapping results so Roger can inspect region under PacBio amplicon primers
	* this region gets trimmed in subsequent steps
	*/

	publishDir "${params.results}/mapped", mode: 'copy'

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
	* A version of this is described in https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/
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
	// Note that this workflow is non-deterministic - the sequence IDs will vary from run-to-run
	// This means that you cannot compare the same data processed through this workflow at different times
	// Instead, reprocess the entire data. Typically this is good as you will be adding new
	// data to an analysis, and the more samples that are included, the greater the likelihood of finding
	// interesting cDNA extensions and novel alleles
	
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
	// This is a major refactoring of the initial allele discovery workflow that avoids 
	// the need for MUSCLE, Clustal Omega, and other tools and is much, much faster
	// TODO : improve handling of adjacent match annotations in CIGAR strings
	
	publishDir params.results, mode: 'copy'

	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(putative_alleles)
	path(ipd_genbank)
	
	output:
	path("${params.experiment_number}-gdna_ipd.gbk"), emit: gdna_ipd
	path("${params.experiment_number}-gdna_extend.gbk"), emit: gdna_extend
	path("${params.experiment_number}-cdna_extend.gbk"), emit: cdna_extend
	path("${params.experiment_number}-novel.gbk"), emit: novel
	path("${params.experiment_number}-mystery.gbk"), emit: mystery
	path("${params.experiment_number}-all.gbk"), emit: all_genbank
	path("${params.experiment_number}-all.fasta"), emit: all_fasta

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
	touch ${params.experiment_number}-gdna_ipd.gbk \
	${params.experiment_number}-gdna_extend.gbk \
	${params.experiment_number}-cdna_extend.gbk \
	${params.experiment_number}-novel.gbk \
	${params.experiment_number}-mystery.gbk \
	${params.experiment_number}-all.gbk \
	${params.experiment_number}-all.fasta

	# classify putative alleles
	classify_putative.py ${putative_alleles} ${params.experiment_number}
	"""

}

process GENOTYPE_PUTATIVE {
	
	// Rename putative allele FASTA files with nearest IPD match
	// Map each animal's clusters to the putative allele FASTA file
	// This shows which animals share a particular putative sequence
	// Filter mappings to only include reads that have an exact match to gDNA sequences (no mismatches).
	// Use the `filterlines.sh` tool to filter SAM files by NM:i:0 tag and `filterbyname.sh` to filter out
	// FASTA sequences that have matches in the filtered SAM files.
	// Then output a CSV version of the genotyping SAM file for use in Excel Pivot Tables.
	
	publishDir params.results, mode: 'copy'

	cpus 40
	errorStrategy 'retry'
	maxRetries 4
	
	input:
		path(clipped_cluster_fasta)
		path(all_putative_fasta)
		
	output:
		path("${params.experiment_number}-genotypes.sam"), emit: genotypes_sam
		path("${params.experiment_number}-genotypes.csv"), emit: csv
	
	script:
	"""
	# replace sequence names in all_fasta with sequence descriptions
	create_renamed_fasta.py ${all_putative_fasta}
	
	# merge renamed source FASTA into single file
	# retain only animal name in sequence headers
	cat *_filtered.fasta | sed 's/_.*//' > spades.fasta

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
		out=${params.experiment_number}-genotypes.sam \
		names=NM:i:0 \
		substring=t \
		include=t

	# create output genotyping CSV file that can be converted to an Excel Pivot Table
	echo "sample,allele" > ${params.experiment_number}-genotypes.csv
	cut -f 1,3 ${params.experiment_number}-genotypes.sam \
	| awk '{OFS="," ; print \$1, \$2}' \
	>> ${params.experiment_number}-genotypes.csv
	"""
	
}

// --------------------------------------------------------------- //


