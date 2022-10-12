#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// input channels
	ch_sample_manifest = Channel
		.fromPath( params.sample_manifest )
		.splitCsv( header: ['bam', 'sample', 'animal'] )
		.map { row -> tuple( "${params.bam_folder}/${row.bam}", row.sample, row.animal ) }
	
	ch_ipd_ref = Channel
		.fromPath( params.ipd_refs )
	
	// Workflow steps
	CONVERT_BAM_TO_FASTQ (
		ch_sample_manifest
	)
	
	ORIENT_FASTQ (
		CONVERT_BAM_TO_FASTQ.out
	)
	
	TRIM_FASTQ (
		ORIENT_FASTQ.out
	)
	
	RUN_PBAA (
		TRIM_FASTQ.out.fastq,
		TRIM_FASTQ.out.index
	)
	
	CLUSTER_PER_SAMPLE (
		RUN_PBAA.out
	)
	
	RENAME_CLUSTERS (
		CLUSTER_PER_SAMPLE.out
	)
	
	MERGE_PER_MAMU_CLUSTERS (
		RENAME_CLUSTERS.out
			.filter { it[2] == "mamu" }
			.map { fasta, sample, animal -> fasta }
			.collect()
	)
	
	MERGE_PER_MAFA_CLUSTERS (
		RENAME_CLUSTERS.out
			.filter { it[2] == "mafa" }
			.map { fasta, sample, animal -> fasta }
			.collect()
	)
	
	SHARED_ANIMALS (
		MERGE_PER_MAMU_CLUSTERS.out
			.mix (
				MERGE_PER_MAFA_CLUSTERS.out
			)
	)
	
	RENAME_PUTATIVE_ALLELE_CLUSTERS (
		SHARED_ANIMALS.out.putative
	)
	
	PARSE_IPD_GENBANK (
		ch_ipd_ref
	)
	
	MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA (
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out,
		PARSE_IPD_GENBANK.out.gdna
	)
	
	FILTER_EXACT_GDNA_MATCHES (
		MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out,
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out
			.map { fasta, animal -> fasta },
		PARSE_IPD_GENBANK.out.gdna
	)
	
	MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches,
		PARSE_IPD_GENBANK.out.cdna
	)
	
	FIND_MUSCLE_CDNA_GDNA_MATCHES (
		MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE.out.merged
	)
	
	RENAME_MUSCLE_CDNA_MATCHES_FASTA (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
			.map { matches, animal -> matches },
		FIND_MUSCLE_CDNA_GDNA_MATCHES.out
	)
	
	EXTRACT_NOVEL_SEQUENCES (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
			.map { matches, animal -> matches },
		RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
	)
	
	MERGE_READS (
		PARSE_IPD_GENBANK.out.gdna,
		EXTRACT_NOVEL_SEQUENCES.out,
		RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
			.map { matches, animal -> matches }
	)
	
	CLUSTAL_ALIGN (
		MERGE_READS.out,
		EXTRACT_NOVEL_SEQUENCES.out
			.map { novel, animal -> novel }
	)
	
	PARSE_DISTANCES (
		CLUSTAL_ALIGN.out.distances,
		EXTRACT_NOVEL_SEQUENCES.out
			.map { novel, animal -> novel },
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
			.map { matches, animal -> matches }
	)
	
	CREATE_GENOTYPING_FASTA (
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out
			.map { clusters, animal -> clusters },
		PARSE_IPD_GENBANK.out.gdna,
		FIND_MUSCLE_CDNA_GDNA_MATCHES.out
			.map { matches, animal -> matches },
		PARSE_DISTANCES.out.closest_matches
	)
	
	GENOTYPE_MAMU_CSS (
		CONVERT_BAM_TO_FASTQ.out,
		CREATE_GENOTYPING_FASTA.out.classified
			.filter { it[1] == 'mamu' }
			.map { fasta, animal -> fasta }
	)
	
	GENOTYPE_MAFA_CSS (
		CONVERT_BAM_TO_FASTQ.out,
		CREATE_GENOTYPING_FASTA.out.classified
			.filter { it[1] == 'mafa' }
			.map { fasta, animal -> fasta }
	)
	
	CREATE_MAMU_GENOTYPING_CSV (
		GENOTYPE_MAMU_CSS.out
			.map { sam, sample, animal -> sam }
			.collect()
	)
	
	CREATE_MAFA_GENOTYPING_CSV (
		GENOTYPE_MAFA_CSS.out
			.map { sam, sample, animal -> sam }
			.collect()
	)
	
	CREATE_GENOTYPING_PIVOT (
		CREATE_MAMU_GENOTYPING_CSV.out
		.mix (
			
			CREATE_MAFA_GENOTYPING_CSV.out
			
		)
	)
	
	PRELIMINARY_EXONERATE_PUTATIVE (
		CREATE_GENOTYPING_FASTA.out.new_allele
	)
	
	PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE (
		PRELIMINARY_EXONERATE_PUTATIVE.out.gff
	)
	
	PRELIMINARY_EXONERATE_MERGE_CDS_PUTATIVE (
		PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE.out.processed_gff
	)
	
}
// --------------------------------------------------------------- //



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Derivative parameters, mostly for making specific results folders
params.raw_fastqs = params.results + "/" + "00-fastq"
params.orient_fastq = params.results + "/" + "01-orient-fastq"
params.trimmed_fastq = params.results + "/" + "02-trim-fastq"
params.pbaa_clusters = params.results + "/" + "03-pbaa/"
params.sample_clusters = params.results + "/" + "04-cluster_per_sample"
params.merged_clusters = params.results + "/" + "05-merged_clusters"
params.shared_clusters = params.results + "/" + "06-shared_clusters"
params.ipd_ref_sep = params.results + "/" + "ipd_ref_separate"
params.gdna_identical = params.results + "/" + "07-gdna-identical"
params.cdna_identical = params.results + "/" + "08-muscle-cdna-identical"
params.novel_alleles = params.results + "/" + "09-novel"
params.genotyping = params.results + "/" + "10-genotyping"
params.genotyping_sams = params.genotyping + "/" + "genotyping_sams"
params.putative_new = params.results + "/" + "11-putative-new-alleles"

params.pbaa_guides = params.pbaa_resources + "/" + "*.fasta"
params.ipd_refs = params.classify_resources + "/" + "*.gbk"
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process CONVERT_BAM_TO_FASTQ {
	
	// Demultiplexed PacBio CCS reads are provided as unmapped BAM files
	// Convert these CCS files to FASTQ format, renaming FASTQ to use sample names
	// The FASTQ filenames need to contain the locus (dpa, dpb, etc.) in the filenames
	// so the appropriate processing parameters can be applied during orienting, trimming, and pbaa
	
	tag "${sample}"
	publishDir params.raw_fastqs, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(bam), val(sample), val(animal)
	
	output:
	tuple path("*.fastq.gz"), val(sample), val(animal)
	
	script:
	basename = bam.getName()
	if( bam.empty() )
		error "The file ${basename} is not present in ${params.bam_folder}"
	else
		"""
		samtools bam2fq ${bam} | gzip > ${sample}.fastq.gz
		"""

}

process ORIENT_FASTQ {
	
	// use vsearch orient command to ensure reads are all in the same orientation
	
	tag "${sample}"
	publishDir params.orient_fastq, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	output:
	tuple path("*.fastq"), val(sample), val(animal)
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		vsearch --orient ${fastq} \
		--db ${params.dpa_reference} \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		vsearch --orient ${fastq} \
		--db ${params.dpb_reference} \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		vsearch --orient ${fastq} \
		--db ${params.dqa_reference} \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		vsearch --orient ${fastq} \
		--db ${params.dqb_reference} \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		vsearch --orient ${fastq} \
		--db ${params.drb_reference} \
		--fastqout ${sample}.fastq
		"""
	
}

process TRIM_FASTQ {
	
	// use bbduk to remove 5' and 3' primer sequences
	// failing to remove these sequences leads to artificial alleles with primer sequences at ends
	// need to trim to minlength or else short sequences choke pbaa
	// since sequences are oriented, only need to trim ends in one orientation
	
	tag "${sample}"
	publishDir params.trimmed_fastq, mode: 'copy'
	
	cpus 1
	memory '2.5 GB'
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	output:
	tuple path("*.fastq"), val(sample), val(animal), emit: fastq
	path("*.fai"), emit: index
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		bbduk.sh -Xmx1g \
		int=f \
		in=${fastq} \
		literal=${params.dpa_forward_primers} \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=${params.dpa_minimum_length} \
		out=stdout.fastq \
		| bbduk.sh -Xmx1g int=f in=stdin.fastq \
		literal=${params.dpa_reverse_primers} \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=${params.dpa_minimum_length} \
		out=${sample}_trimmed.fastq && \
		samtools fqidx ${sample}_trimmed.fastq
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		bbduk.sh -Xmx1g \
		int=f \
		in=${fastq} \
		literal=${params.dpb_forward_primers} \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=${params.dpb_minimum_length} \
		out=stdout.fastq \
		| bbduk.sh -Xmx1g int=f in=stdin.fastq \
		literal=${params.dpb_reverse_primers} \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=${params.dpb_minimum_length} \
		out=${sample}_trimmed.fastq && \
		samtools fqidx ${sample}_trimmed.fastq
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		bbduk.sh -Xmx1g \
		int=f \
		in=${fastq} \
		literal=${params.dqa_forward_primers} \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=${params.dqa_minimum_length} \
		out=stdout.fastq \
		| bbduk.sh -Xmx1g int=f in=stdin.fastq \
		literal=${params.dqa_reverse_primers} \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=${params.dqa_minimum_length} \
		out=${sample}_trimmed.fastq && \
		samtools fqidx ${sample}_trimmed.fastq
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		bbduk.sh -Xmx1g \
		int=f \
		in=${fastq} \
		literal=${params.dqb_forward_primers} \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=${params.dqb_minimum_length} \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=${params.dqb_reverse_primers} \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=${params.dqb_minimum_length} \
		out=${sample}_trimmed.fastq && \
		samtools fqidx ${sample}_trimmed.fastq
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		bbduk.sh -Xmx1g \
		int=f \
		in=${fastq} \
		literal=${params.drb_forward_primers} \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=${params.drb_minimum_length} \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=${params.drb_reverse_primers} \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=${params.drb_minimum_length} \
		out=${sample}_trimmed.fastq && \
		samtools fqidx ${sample}_trimmed.fastq
		"""
}

process RUN_PBAA {
	
	// run pbaa clustering on each FASTQ that has been trimmed
	// suppress stderr so you can see which samples fail
	
	tag "${sample}"
	publishDir params.pbaa_clusters, pattern: "*.fasta", mode: 'copy'
	
	cpus 2
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	path(index)
	
	output:
	tuple path("*.fasta"), val(sample), val(animal)
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		samtools faidx ${params.dpa_guide_fasta} && \
		/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		${params.dpa_guide_fasta} \
		${fastq} \
		"${sample}" 2> /dev/null && \
		mv ${sample}_passed_cluster_sequences.fasta ${sample}.fasta
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		samtools faidx ${params.dpb_guide_fasta} && \
		/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		${params.dpb_guide_fasta} \
		${fastq} \
		"${sample}" 2> /dev/null && \
		mv ${sample}_passed_cluster_sequences.fasta ${sample}.fasta
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		samtools faidx ${params.dqa_guide_fasta} && \
		/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		${params.dqa_guide_fasta} \
		${fastq} \
		"${sample}" 2> /dev/null && \
		mv ${sample}_passed_cluster_sequences.fasta ${sample}.fasta
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		samtools faidx ${params.dqb_guide_fasta} && \
		/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		${params.dqb_guide_fasta} \
		${fastq} \
		"${sample}" 2> /dev/null && \
		mv ${sample}_passed_cluster_sequences.fasta ${sample}.fasta
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		samtools faidx ${params.drb_guide_fasta} && \
		/miniconda2/bin/pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		${params.drb_guide_fasta} \
		${fastq} \
		"${sample}" 2> /dev/null && \
		mv ${sample}_passed_cluster_sequences.fasta ${sample}.fasta
		"""

}

process CLUSTER_PER_SAMPLE {
	
	// run bbbmap dedupe.sh on each sample, finding exact matches and absorbing containments.
	// use PacBio clustering parameters from https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/dedupe-guide/ but change edit distance to 0.
	// This should restrict clustering to identical sequences. 
	// outbest parameter is not documented in dedupe.sh but is used in the example from the PacBio clustering and saves one representative sequence per cluster.
	// Since the sequences in the cluster are identical and have an edit distance of 0, this exemplar sequence should be sufficient.
	// 
	// Some samples might not have output. Create empty output file and overwrite if data exists.
	
	tag "${sample}"
	publishDir params.sample_clusters, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(sample), val(animal)
	
	output:
	tuple path("*_clustered.fasta.gz"), val(sample), val(animal)
	
	script:
	"""
	gzip < /dev/null > ${sample}.fasta.gz
	dedupe.sh -Xmx1g ow=t \
	in=${sample}.fasta \
	outbest=${sample}_clustered.fasta.gz \
	fo c \
	threads=${task.cpus}
	"""

}

process RENAME_CLUSTERS {
	
	// prepend sample name to cluster names. This makes it easier to track which sequences are associated with which clusters.
	
	tag "${sample}"
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(sample), val(animal)
	
	output:
	tuple path("*.fasta.gz"), val(sample), val(animal)
	
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

process MERGE_PER_MAMU_CLUSTERS {
	
	// gymnastics with per-sample merged FASTA to merge mamu samples after all have been created - need to wait until rename_clusters completes before running this step
	
	publishDir params.merged_clusters, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(mamu_list)
	
	output:
	tuple path("*fasta.gz"), val("mamu")
	
	script:
	"""
	zcat *.fasta.gz | gzip > merged_mamu_clusters.fasta.gz
	"""

}

process MERGE_PER_MAFA_CLUSTERS {
	
	// gymnastics with per-sample merged FASTA to merge mafa samples after all have been created - need to wait until rename_clusters completes before running this step
	
	publishDir params.merged_clusters, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(mafa_list)
	
	output:
	tuple path("*fasta.gz"), val("mafa")
	
	script:
	"""
	zcat *.fasta.gz | gzip > merged_mafa_clusters.fasta.gz
	"""

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
	
	tag "${animal}"
	publishDir params.shared_clusters, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(animal)
	
	output:
	tuple path("*shared.fasta"), val(animal), emit: shared
	tuple path("*all.fasta"), val(animal), emit: all
	tuple path("*unique.fasta"), val(animal), emit: unique
	tuple path("*putative_alleles_temp.fasta"), val(animal), emit: putative
	
	
	script:
	"""
	dedupe.sh -Xmx1g \
	in=${fasta} \
	outbest=${animal}_all.fasta \
	am=t ac=f arc=t fo c fcc nam=4 threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=${fasta} \
	out=${animal}_unique.fasta \
	am=t ac=f arc=t fo fcc uniqueonly=t threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=${animal}_all.fasta,${animal}_unique.fasta \
	out=${animal}_shared.fasta \
	ac=f uniqueonly=t threads=${task.cpus}
	
	dedupe.sh -Xmx1g \
	in=${fasta} \
	out=${animal}_putative_alleles_temp.fasta \
	ac=t threads=${task.cpus}
	"""

}

process RENAME_PUTATIVE_ALLELE_CLUSTERS {
	
	// add integer FASTA id to gdna_match FASTA name
	// this simplifies downstream analyses, where the FASTA description can be complex and the FASTA id is simple
	// since previous steps are nondeterminisitic, add timestamp to id too
	// this makes it easier to clarify which files are used for genotyping
	
	tag "${animal}"
	publishDir params.shared_clusters, mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(putative), val(animal)
	
	output:
	tuple path("*putative_alleles.fasta"), val(animal)
	
	script:
	"""
	#!/usr/bin/env python3
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	from datetime import datetime
	
	with open("${animal}" + "_putative_alleles.fasta", "w") as handle:
		for idx, record in enumerate(SeqIO.parse("${putative}", "fasta")):
			now = datetime.now()
			time = now.strftime("%Y%m%d%H%M%S")
			record.id = str(time) + '-' + str(idx)
			SeqIO.write(record, handle, "fasta")
	"""

}

process PARSE_IPD_GENBANK {
	
	// use biopython to take IPD Genbank file and separate into FASTA
	// of genomic DNA (intron containing) sequences and cDNA sequences (no intron)
	
	tag "${animal_name}"
	publishDir params.ipd_ref_sep, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(guide_fasta)
	
	output:
	path("*gdna_reference.fasta"), emit: gdna
	path("*cdna_reference.fasta"), emit: cdna
	
	script:
	animal_name = guide_fasta.simpleName.substring(8,12)
	
	"""
	ipd_to_gdna_cdna.py ${guide_fasta} ${animal_name}
	"""

}

process MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA {
	
	// identify putative alleles whose sequences match full-length gDNA sequences already in IPD.
	// save SAM file of reads that map to known IPD alleles and FASTA.gz file of reads that do not.
	
	tag "${putative_animal}"
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	putative_animal == gdna.simpleName.substring(0,4)
	
	input:
	tuple path(putative), val(putative_animal)
	each path(gdna) 
	
	output:
	tuple path("*all_mappings.sam"), val(putative_animal)
	
	script:
	"""
	minimap2 \
	-t ${task.cpus} \
	${gdna} \
	${putative} \
	-ax splice -N 10000 \
	> ${putative_animal}_all_mappings.sam
	"""
	

}

process FILTER_EXACT_GDNA_MATCHES {
	
	// filter mappings to only those that have NM:i:0 (no mismatches)
	// use filterlines.sh tool
	
	tag "${putative_animal}"
	publishDir params.gdna_identical, pattern: '*gdna_match.sam', mode: 'copy'
	publishDir params.cdna_identical, pattern: '*no-gdna_match.fasta', mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	putative_animal == fasta.simpleName.substring(0,4) && putative_animal == gdna.simpleName.substring(0,4)
	
	input:
	tuple path(sam), val(putative_animal)
	each path(fasta)
	each path(gdna)
	
	output:
	tuple path("*_gdna_match.sam"), val(putative_animal), emit: gdna_matches
	tuple path("*_no-gdna_match.fasta"), val(putative_animal), emit: cdna_matches
	
	script:
	"""
	filterlines.sh \
	in=${sam} \
	out=${putative_animal}_gdna_match.sam \
	names=NM:i:0 \
	substring=t \
	include=t
	
	filterbyname.sh \
	in=${fasta} \
	names=${putative_animal}_gdna_match.sam \
	out=${putative_animal}_no-gdna_match.fasta
	"""

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
	
	when:
	novel_animal == cdna_matches.simpleName.substring(0,4) && novel_animal == gdna_ref.simpleName.substring(0,4)
	
	input:
	each path(gdna_ref)
	tuple path(novel), val(novel_animal)
	each path(cdna_matches)
	
	output:
	tuple path("*reads.fasta"), val(novel_animal)
	
	script:
	"""
	cat ${gdna_ref} ${novel} ${cdna_matches} > ${novel_animal}_reads.fasta
	"""

}

process CLUSTAL_ALIGN {
	
	// align reads with clustal omega
	// generate alignment in fasta format and distance matrix
	// distance matrix can be used to parse closest matches to known alleles
	
	tag "${reads_animal}"
	publishDir params.novel_alleles, pattern: '*distances.txt', mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	reads_animal == novel.simpleName.substring(0,4)
	
	input:
	tuple path(reads), val(reads_animal)
	each path(novel)
	
	output:
	tuple path("*aligned.fasta"), val(reads_animal)
	tuple path("*distances.txt"), val(reads_animal), emit: distances
	
	script:
	"""
	clustalo \
	--infile=${reads} \
	--outfile=${reads_animal}_aligned.fasta \
	--distmat-out=${reads_animal}_distances.txt \
	--threads=${task.cpus} \
	--full
	"""

}

process PARSE_DISTANCES {
	
	// parse clustal omega distances to find alleles in reference database and cDNA extensions that are nearest neighbors to novel sequences
	
	tag "${clustal_animal}"
	
	errorStrategy 'retry'
	maxRetries 4
	
	when:
	clustal_animal == novel.simpleName.substring(0,4) && clustal_animal == cdna_matches.simpleName.substring(0,4)
	
	input:
	tuple path(distances), val(clustal_animal)
	each path(novel)
	each path(cdna_matches)
	
	output:
	tuple path("*novel_closest_matches.xlsx"), val(clustal_animal), emit: closest_matches
	tuple path("*distances_tmp.txt"), val(clustal_animal)
	
	script:
	"""
	parse_clustalo_distances.py ${clustal_animal} ${novel} ${cdna_matches}
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

process GENOTYPE_MAMU_CSS {
	
	// identify reference sequences that are fully contained in CCS reads
	// use trimmed full-length reads for mapping
	// use higher minratio and smaller maxindel per emails with Brian Bushnell
	
	tag "${sample}"
	publishDir params.genotyping_sams, pattern: '*.sam', mode: 'copy'
	
	cpus 4
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	each path(classified)
	
	output:
	tuple path('*.sam'), val(sample), val(animal)
	
	script:
	"""
	minimap2 \
	${fastq} \
	${classified} \
	-ax map-hifi --eqx -t 3 \
	| reformat.sh \
	in=stdin.sam \
	out=${sample}.sam \
	ref=${classified} \
	noheader=t \
	editfilter=0 \
	threads=1 \
	ow=t
	"""

}

process GENOTYPE_MAFA_CSS {
	
	// identify reference sequences that are fully contained in CCS reads
	// use trimmed full-length reads for mapping
	// use higher minratio and smaller maxindel per emails with Brian Bushnell
	
	tag "${sample}"
	publishDir params.genotyping_sams, pattern: '*.sam', mode: 'copy'
	
	cpus 4
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	each path(classified)
	
	output:
	tuple path('*.sam'), val(sample), val(animal)
	
	script:
	"""
	minimap2 \
	${fastq} \
	${classified} \
	-ax map-hifi --eqx -t 3 \
	| reformat.sh \
	in=stdin.sam \
	out=${sample}.sam \
	ref=${classified} \
	noheader=t \
	editfilter=0 \
	threads=1 \
	ow=t
	"""

}

process CREATE_MAMU_GENOTYPING_CSV {
	
	// create CSV file containing sample name, allele name, and count of reads matching allele
	
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	path(sam_list)
	
	output:
	tuple path("*.csv"), val("mamu")
	
	script:
	sam_string = sam_list.join(",")
	
	"""
	create_genotyping_csv.py "mamu" ${sam_string}
	"""

}


process CREATE_MAFA_GENOTYPING_CSV {
	
	// create CSV file containing sample name, allele name, and count of reads matching allele

	errorStrategy 'retry'
	maxRetries 4
		
	input:
	path(sam_list)
	
	output:
	tuple path("*.csv"), val("mafa")
	
	script:
	sam_string = sam_list.join(",")
	
	"""
	create_genotyping_csv.py "mafa" ${sam_string}
	"""

}


process CREATE_GENOTYPING_PIVOT {
	
	// make Excel-formatted pivot table from genotyping CSV
	
	publishDir params.genotyping, pattern: '*.sam', mode: 'copy'
	
	errorStrategy 'retry'
	maxRetries 4
	
	tag "${animal}"
	publishDir params.genotyping, mode: 'move'
	
	input:
	tuple path(csv), val(animal)
	
	output:
	
	script:
	"""
	genotyping.py ${animal}
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


