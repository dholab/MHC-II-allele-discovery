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
	
	ch_ipd_ref = Channel
		.fromPath( params.ipd_refs )
	
	// Workflow steps
	
	ORIENT_FASTQ (
		ch_sample_manifest
	)

	MAP_FASTQ (
		ORIENT_FASTQ.out.oriented_fasta
	)

	TRIM_TO_PACBIO_AMPLICONS (
		MAP_FASTQ.out.mapped_bam
	)

	FILTER_HARD_CLIPPED_AMPLICONS (
		TRIM_TO_PACBIO_AMPLICONS.out.trimmed_bam
	)
	
	RENAME_CLUSTERS (
		FILTER_HARD_CLIPPED_AMPLICONS.out.clipped_fasta
	)
	
	MERGE_PER_SAMPLE_CLUSTERS (
		RENAME_CLUSTERS.out.renamed_fasta
			.map { fasta, sample-> fasta }
			.collect()
	)
	
	SHARED_ANIMALS (
		MERGE_PER_SAMPLE_CLUSTERS.out.merged_fasta
	)
	
	RENAME_PUTATIVE_ALLELE_CLUSTERS (
		SHARED_ANIMALS.out.putative
	)
	
	PARSE_IPD_GENBANK (
		ch_ipd_ref
	)
	
	MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA (
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out.renamed_clusters,
		PARSE_IPD_GENBANK.out.ipd_gdna
	)
	
	FILTER_EXACT_GDNA_MATCHES (
		MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out.all_mappings
	)
	
	DEFINE_CDNA_MATCHES_AND_NOVELS (
		FILTER_EXACT_GDNA_MATCHES.out.no_gdna_matches_fasta,
		PARSE_IPD_GENBANK.out.ipd_cdna
	)

	// MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE (
	// 	FILTER_EXACT_GDNA_MATCHES.out.cdna_matches,
	// 	PARSE_IPD_GENBANK.out.cdna
	// )
	
	// FIND_MUSCLE_CDNA_GDNA_MATCHES (
	// 	MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE.out.merged
	// )
	
	// RENAME_MUSCLE_CDNA_MATCHES_FASTA (
	// 	FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
	// 		.map { matches, animal -> matches },
	// 	FIND_MUSCLE_CDNA_GDNA_MATCHES.out
	// )
	
	// EXTRACT_NOVEL_SEQUENCES (
	// 	FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
	// 		.map { matches, animal -> matches },
	// 	RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
	// )
	
	MERGE_READS (
		PARSE_IPD_GENBANK.out.ipd_gdna,
		DEFINE_CDNA_MATCHES_AND_NOVELS.out.cdna,
		DEFINE_CDNA_MATCHES_AND_NOVELS.out.novel
	)
	
	CLUSTAL_ALIGN (
		MERGE_READS.out
	)
	
	PARSE_DISTANCES (
		CLUSTAL_ALIGN.out.distances,
		DEFINE_CDNA_MATCHES_AND_NOVELS.out.novel,
		DEFINE_CDNA_MATCHES_AND_NOVELS.out.cdna
	)
	
	// CREATE_GENOTYPING_FASTA (
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out
	// 		.map { clusters, animal -> clusters },
	// 	PARSE_IPD_GENBANK.out.gdna,
	// 	FIND_MUSCLE_CDNA_GDNA_MATCHES.out
	// 		.map { matches, animal -> matches },
	// 	PARSE_DISTANCES.out.closest_matches
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
params.orient_fastq = params.results + "/" + "00-orient-fastq"
params.map_fastq = params.results + "/" + "01-map-fastq"
params.trim_to_pacbio_amplicons = params.results + "/" + "02-trim-to-pacbio-amplicons"
params.filter_hard_clipped_amplicons = params.results + "/" + "03-filter-hard-clipped-amplicons"
params.sample_clusters = params.results + "/" + "04-cluster_per_sample"
params.merged_clusters = params.results + "/" + "05-merged_clusters"
params.shared_clusters = params.results + "/" + "06-shared_clusters"
params.ipd_ref_sep = params.results + "/" + "ipd_ref_separate"
params.gdna_identical = params.results + "/" + "07-gdna-identical"
params.cdna_identical = params.results + "/" + "08-muscle-cdna-identical"
params.novel_alleles = params.results + "/" + "09-novel"
params.genotyping = params.results + "/" + "10-genotyping"
params.putative_new = params.results + "/" + "11-putative-new-alleles"
params.ipd_refs = params.classify_resources + "/" + "*.gbk"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process ORIENT_FASTQ {
	/*
	* Process: ORIENT_FASTQ
	* Description: Use vsearch orient command to ensure reads are all in the same orientation.
	* This avoids complex reverse complementing some contigs in subsequent steps.
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	publishDir params.orient_fastq, mode: 'copy'
	
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

process MAP_FASTQ {
	
	/*
	* Process: MAP_FASTQ
	* Description: Use minimap2 to map SPAdes contigs to combined DPA, DPB, DQA, DQB, DRB reference FASTA.
	* This reference sequence contains 1kb of flanking sequence at end of each gene,
	* which is needed to ensure that indels near the PacBio amplicon ends map correctly.
	* After mapping, use samtools to sort the minimap2 output and convert to BAM file.
	*/

	// Set process label and publish directory for output files
	tag "${sample}"
	publishDir params.map_fastq, mode: 'copy'
	
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
	publishDir params.trim_to_pacbio_amplicons, mode: 'copy'
	
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
	 *   - publishDir: directory to publish output files to
	 */
	tag "${sample}"
	publishDir params.filter_hard_clipped_amplicons, mode: 'copy'

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

	publishDir params.merged_clusters, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(mamu_list)
	
	output:
	path("merged_clusters.fasta.gz"), emit: merged_fasta
	
	script:
	"""
	zcat *.fasta.gz | gzip > merged_clusters.fasta.gz
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

	publishDir params.shared_clusters, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta)
	
	output:
	tuple path("putative_alleles_temp.fasta"), emit: putative
	
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
	in=${fasta} \
	out=putative_alleles_temp.fasta \
	ac=t threads=${task.cpus}
	"""
	
}

process RENAME_PUTATIVE_ALLELE_CLUSTERS {
	
	// Add integer FASTA id to gdna_match FASTA name to simplify downstream analyses.
	// A timestamp is also added to the id to help identify the files used for genotyping.
	
	publishDir params.shared_clusters, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(putative)
	
	output:
	tuple path("putative_alleles.fasta"), emit: renamed_clusters
	
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

process PARSE_IPD_GENBANK {
	
	// Use Biopython to extract genomic DNA and cDNA sequences from an IPD Genbank file. The output files are saved as FASTA format..
	
	publishDir params.ipd_ref_sep, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(guide_fasta)
	
	output:
	path("gdna_reference.fasta"), emit: ipd_gdna
	path("cdna_reference.fasta"), emit: ipd_cdna
	
	"""
	ipd_to_gdna_cdna.py ${guide_fasta} ${animal_name}
	"""

}

process MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA {
	
	// Use minimap2 to map putative alleles whose sequences match full-length genomic DNA (gDNA) sequences already in the IPD database.
	// Save a SAM file of reads that map to known IPD alleles and a gzipped FASTA file of reads that do not.
	
	cpus 3
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(putative_alleles)
	path(ipd_gdna) 
	
	output:
	path("all_mappings.sam"), emit: all_mappings
	
	script:
	"""
	minimap2 -t ${task.cpus} \
	${ipd_gdna} \
	${putative_alleles} \
	-ax splice -N 10000 \
	> all_mappings.sam
	"""
}

process FILTER_EXACT_GDNA_MATCHES {
	
	// Filter mappings to only include reads that have an exact match to gDNA sequences (no mismatches).
	// Use the `filterlines.sh` tool to filter SAM files by NM:i:0 tag and `filterbyname.sh` to filter out
	// FASTA sequences that have matches in the filtered SAM files.
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
		path(all_mappings_sam)
	
	output:
		path("gdna_match.sam"), emit: gdna_matches_sam
		path("*_no-gdna_match.fasta"), emit: no_gdna_matches_fasta
	
	script:
	"""
	filterlines.sh \
		in=${all_mappings_sam} \
		out=gdna_match.sam \
		names=NM:i:0 \
		substring=t \
		include=t
	
	filterbyname.sh \
		in=${fasta} \
		names=gdna_match.sam \
		out=no-gdna_match.fasta
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
	
	tag "${putative_animal}"
	publishDir params.cdna_identical, pattern: '*merged.aln', mode: 'copy'
	publishDir params.cdna_identical, pattern: '*gdna_single_temp.fasta', mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	putative_animal == cdna.simpleName.substring(0,4)
	
	input:
	tuple path(cdna_matches), val(putative_animal)
	each path(cdna)
	
	output:
	tuple path("*merged.aln"), val(putative_animal), emit: merged
	tuple path("*gdna_single_temp.fasta"), val(putative_animal)
	
	script:
	"""
	muscle_cdna_mapping.py ${cdna_matches} ${putative_animal} ${cdna}
	"""
	
}

process FIND_MUSCLE_CDNA_GDNA_MATCHES {
	
	// run awk to find gdna sequences that fully match sequence of cdna
	
	tag "${putative_animal}"
	publishDir params.cdna_identical, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(merged), val(putative_animal)
	
	output:
	tuple path("*matches.aln"), val(putative_animal)
	
	shell:
	'''
	awk \'{{if( $4 == $5  ) print $0}}\' !{merged} > !{putative_animal}_matches.aln
	'''

}

process RENAME_MUSCLE_CDNA_MATCHES_FASTA {
	
	// input gDNA sequences
	// change the name of the sequence header
	// write output sequences that match cDNAs to new file
	
	tag "${gdna_animal}"
	publishDir params.cdna_identical, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	gdna_animal == cdna_matches.simpleName.substring(0,4)
	
	input:
	each path(cdna_matches)
	tuple path(matches_aln), val(gdna_animal)
	
	output:
	tuple path("*cdna_matches.fasta"), val(cdna_animal)
	
	script:
	cdna_animal = cdna_matches.simpleName.substring(0,4)
	"""
	#!/usr/bin/env python3
	
	import subprocess
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	import csv
	
	# create output file in case it is empty
	subprocess.run('touch ${cdna_animal}_cdna_matches.fasta', shell=True)
	
	# read input FASTA line-by-line
	for record in SeqIO.parse("${cdna_matches}", "fasta"):
	
		# parse file with gdna sequences that match cdna sequences
		with open("${matches_aln}") as tsvfile:
			reader = csv.reader(tsvfile, delimiter='\t')
			for row in reader:
	
				# test if name of sequence in cdna match file matches gdna sequence name
				if row[0] == record.name:
					# update name of sequence in output file
					record.description = row[1] + '|' + row[0]
	
					# write to file
					with open("${cdna_animal}_cdna_matches.fasta", "a") as handle:
						SeqIO.write(record, handle, "fasta")
	"""
	
}

process EXTRACT_NOVEL_SEQUENCES {
	
	// the previous steps processed clusters that perfectly match known IPD cDNA sequences.
	// Because the sequences in this experiment are supported by PacBio cDNA and pbaa gDNA sequences,
	// we can submit them to IPD even if they don't match a current IPD cDNA sequence.
	// 
	// Map gDNA sequences to cDNAs and report reads that don't map
	
	tag "${cdna_animal}"
	publishDir params.novel_alleles, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	no_gdna_match.simpleName.substring(0,4) == cdna_animal
	
	input:
	each path(no_gdna_match)
	tuple path(cdna_matches), val(cdna_animal)
	
	output:
	tuple path("*novel.fasta"), val(cdna_animal)
	
	script:
	no_gdna_animal = no_gdna_match.simpleName.substring(0,4)
	
	if( file(cdna_matches).isEmpty() )
		"""
		touch ${no_gdna_animal}_novel.fasta
		"""
	else
		"""
		mapPacBio.sh in=${no_gdna_match} ref=${cdna_matches} outu=${no_gdna_animal}_novel.fasta subfilter=0
		"""

}

process MERGE_READS {
	
	// concatenate files to align with clustal omega
	
	tag "${novel_animal}"
	
	errorStrategy 'retry'
	maxRetries 4
	input:
	
	path(gdna_ref)
	path(cdna_match)
	path(novel)
	
	output:
	path("*reads.fasta")
	
	script:
	"""
	cat ${gdna_ref} ${cdna_match} ${novel} > merged_pre_clustal_reads.fasta
	"""

}

process CLUSTAL_ALIGN {
	
	// align reads with clustal omega
	// generate alignment in fasta format and distance matrix
	// distance matrix can be used to parse closest matches to known alleles
	cpus 40
  	executor 'local'

	tag "${reads_animal}"
	publishDir params.novel_alleles, pattern: '*distances.txt', mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(merged)
	
	output:
	path("*aligned.fasta"), emit: aligned
	path("*distances.txt"), emit: distances

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

process PARSE_DISTANCES {
	
	// parse clustal omega distances to find alleles in reference database and cDNA extensions that are nearest neighbors to novel sequences
	
	tag "${clustal_animal}"
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(distances)
	path(novel)
	path(cdna_matches)
	
	output:
	path("*novel_closest_matches.xlsx"), emit: closest_matches
	path("*distances_tmp.txt"), emit: distances
	
	script:
	"""
	parse_clustalo_distances.py ${novel} ${cdna_matches}
	"""

}

process CREATE_GENOTYPING_FASTA {
	
	// to check accuracy of novel sequences, want to genotype IPD gDNA matches, cDNA matches, and novel sequences
	// need to create FASTA file that contains all of these putative_alleles along with their classification
	// 
	// 27319 - create FASTA file that only has cDNA extensions and novel alleles
	
	tag "${closest_animal}"
	publishDir params.genotyping, pattern: '*_classified.fasta', mode: 'copy'
	publishDir params.putative_new, pattern: '*_putative.fasta', mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	putative.simpleName.substring(0,4) == gdna_ref.simpleName.substring(0,4) && putative.simpleName.substring(0,4) == matches.simpleName.substring(0,4) && putative.simpleName.substring(0,4) == closest_animal
	
	input:
	each path(putative)
	each path(gdna_ref)
	each path(matches)
	tuple path(closest_matches), val(closest_animal)
	
	output:
	tuple path("*_classified.fasta"), val(closest_animal), emit: classified
	tuple path("*_putative.fasta"), val(closest_animal), emit: new_allele
	
	script:
	"""
	create_genotyping_fasta.py ${closest_animal} ${putative} ${gdna_ref} ${matches}
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
	
	tag "${animal}"
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
	
	tag "${animal}"
	
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
	
	tag "${animal}"
	publishDir params.putative_new, mode: 'move'
	
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


