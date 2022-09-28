#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


// Specifying additional input file types
params.pbaa_guides = params.pbaa_resources + "/" + "*.fasta"
params.ipd_refs = params.classify_resources + "/" + "*.gbk"


workflow {
	
	// input channels
	ch_sample_manifest = Channel
		.fromPath( params.sample_manifest )
		.splitCsv( sep: \t, header: ['bam', 'sample', 'animal'] )
		.map { row -> tuple( 
			file("${params.bam_folder}/${row.bam}"), row.sample, row.animal 
		)}
	
	ch_pbaa_guides = Channel
		.fromPath( params.pbaa_guides )
	
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
	
	INDEX_FASTQ (
		TRIM_FASTQ.out
	)
	
	INDEX_GUIDE_FASTA (
		ch_pbaa_guides
	)
	
	RUN_PBAA (
		INDEX_FASTQ.out.fastq,
		INDEX_FASTQ.out.index,
		ch_pbaa_guides.collect(),
		INDEX_GUIDE_FASTA.out.collect()
	)
	
	RENAME_FOR_LAA (
		RUN_PBAA.out
	)
	
	CLUSTER_PER_SAMPLE (
		RENAME_FOR_LAA.out
	)
	
	RENAME_CLUSTERS (
		CLUSTER_PER_SAMPLE.out
	)
	
	MERGE_PER_MAMU_CLUSTERS (
		RENAME_CLUSTERS.out
		.map { it -> tuple(it.fasta, it.sample, it.animal)}
		.filter { it[2] == "mamu" }
		.map { it[0], it[1], it[2] -> it[0] }
		.collect()
	)
	
	MERGE_PER_MAFA_CLUSTERS (
		RENAME_CLUSTERS.out
		.map { it -> tuple(it.fasta, it.sample, it.animal)}
		.filter { it[2] == "mafa" }
		.map { it[0], it[1], it[2] -> it[0] }
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
		PARSE_IPD_GENBANK.gdna,
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out
	)
	
	FILTER_EXACT_GDNA_MATCHES (
		MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out,
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out,
		PARSE_IPD_GENBANK.out.gdna
	)
	
	MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_match,
		PARSE_IPD_GENBANK.out.cdna
	)
	
	FIND_MUSCLE_CDNA_GDNA_MATCHES (
		MAP_SHARED_CLUSTERS_TO_CDNA_WITH_MUSCLE.out.merged
	)
	
	RENAME_MUSCLE_CDNA_MATCHES_FASTA (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches,
		FIND_MUSCLE_CDNA_GDNA_MATCHES.out
	)
	
	EXTRACT_NOVEL_SEQUENCES (
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches,
		RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
	)
	
	MERGE_READS (
		PARSE_IPD_GENBANK.out.gdna,
		EXTRACT_NOVEL_SEQUENCES.out,
		RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
	)
	
	CLUSTAL_ALIGN (
		MERGE_READS.out,
		EXTRACT_NOVEL_SEQUENCES.out
	)
	
	PARSE_DISTANCES (
		CLUSTAL_ALIGN.out.distances,
		EXTRACT_NOVEL_SEQUENCES.out,
		FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
	)
	
	CREATE_GENOTYPING_FASTA (
		RENAME_PUTATIVE_ALLELE_CLUSTERS.out,
		PARSE_IPD_GENBANK.out.gdna,
		FIND_MUSCLE_CDNA_GDNA_MATCHES.out,
		PARSE_DISTANCES.out.closest_matches
	)
	
	GENOTYPE_CSS (
		CONVERT_BAM_TO_FASTQ.out,
		CREATE_GENOTYPING_FASTA.out.classified
	)
	
	CREATE_MAMU_GENOTYPING_CSV (
		GENOTYPE_CSS.out
		.map { it -> tuple(it.sam, it.sample, it.animal)}
		.filter { it[2] == "mamu" }
		map { it[0], it[1], it[2] -> it[0] }
		.collect()
	)
	
	CREATE_MAFA_GENOTYPING_CSV (
		GENOTYPE_CSS.out
		.map { it -> tuple(it.sam, it.sample, it.animal)}
		.filter { it[2] == "mafa" }
		map { it[0], it[1], it[2] -> it[0] }
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


// Derivative parameters, mostly for making specific results folders
params.raw_fastqs = params.results + "/" + "00-fastq"
params.orient_fastq = params.results + "/" + "01-orient-fastq"
params.trimmed_fastq = params.results + "/" + "02-trim-fastq"
params.pbaa_clusters = params.results + "/" + "04-pbaa/"
params.sample_clusters = params.results + "/" + "05-cluster_per_sample"
params.merged_clusters = params.results + "/" + "06-merged_clusters"
params.shared_clusters = params.results + "/" + "07-merged_clusters"
params.ipd_ref_sep = params.results + "/" + "ipd_ref_separate"
params.gdna_identical = params.results + "/" + "08-gdna-identical"
params.cdna_identical = params.results + "/" + "09-muscle-cdna-identical"
params.novel_alleles = params.results + "/" + "10-novel"
params.genotyping = params.results + "/" + "11-genotyping"
params.putative_new = params.results + "/" + "12-putative-new-alleles"

	
process CONVERT_BAM_TO_FASTQ {
	
	// Demultiplexed PacBio CCS reads are provided as unmapped BAM files
	// Convert these CCS files to FASTQ format, renaming FASTQ to use sample names
	// The FASTQ filenames need to contain the locus (dpa, dpb, etc.) in the filenames
	// so the appropriate processing parameters can be applied during orienting, trimming, and pbaa
	
	tag "${sample}"
	publishDir params.raw_fastqs, mode: 'copy'
	
	cpus 1
	
	input:
	tuple path(bam), val(sample), val(animal)
	
	output:
	tuple path("*.fastq.gz"), val(sample), val(animal)
	
	script:
	"""
	samtools bam2fq ${bam} | gzip > ${sample}.fastq.gz
	"""

}

process ORIENT_FASTQ {
	
	// use vsearch orient command to ensure reads are all in the same orientation
	
	tag "${sample}"
	publishDir params.orient_fastq, mode: 'symlink'
	
	cpus 1
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	output:
	tuple path("*.fastq"), val(sample), val(animal)
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		vsearch --orient ${fastq} \
		--db params.orient.dpa.orient_reference \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		vsearch --orient ${fastq} \
		--db params.orient.dpb.orient_reference \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		vsearch --orient ${fastq} \
		--db params.orient.dqa.orient_reference \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		vsearch --orient ${fastq} \
		--db params.orient.dqb.orient_reference \
		--fastqout ${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		vsearch --orient ${fastq} \
		--db params.orient.drb.orient_reference \
		--fastqout ${sample}.fastq
		"""
	
}

process TRIM_FASTQ {
	
	// use bbduk to remove 5' and 3' primer sequences
	// failing to remove these sequences leads to artificial alleles with primer sequences at ends
	// need to trim to minlength or else short sequences choke pbaa
	// since sequences are oriented, only need to trim ends in one orientation
	
	tag "${sample}"
	publishDir params.trimmed_fastq, mode: 'symlink'
	
	cpus 1
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		bbduk.sh int=f \
		in=${fastq} \
		literal=params.trim.dpa.forward_primers \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=params.trim.dpa.minimum_length \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=params.trim.dpa.reverse_primers \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=params.trim.dpa.minimum_length \
		out=${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		bbduk.sh int=f \
		in=${fastq} \
		literal=params.trim.dpb.forward_primers \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=params.trim.dpb.minimum_length \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=params.trim.dpb.reverse_primers \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=params.trim.dpb.minimum_length \
		out=${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		bbduk.sh int=f \
		in=${fastq} \
		literal=params.trim.dqa.forward_primers \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=params.trim.dqa.minimum_length \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=params.trim.dqa.reverse_primers \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=params.trim.dqa.minimum_length \
		out=${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		bbduk.sh int=f \
		in=${fastq} \
		literal=params.trim.dqb.forward_primers \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=params.trim.dqb.minimum_length \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=params.trim.dqb.reverse_primers \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=params.trim.dqb.minimum_length \
		out=${sample}.fastq
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		bbduk.sh int=f \
		in=${fastq} \
		literal=params.trim.drb.forward_primers \
		restrictleft=50 ktrim=l k=8 qin=33 \
		minlength=params.trim.drb.minimum_length \
		out=stdout.fastq \
		| bbduk.sh int=f in=stdin.fastq \
		literal=params.trim.drb.reverse_primers \
		restrictright=50 ktrim=r k=8 qin=33 \
		minlength=params.trim.drb.minimum_length \
		out=${sample}.fastq
		"""
	
	

}

process INDEX_FASTQ {
	
	// make index from each FASTQ file
	// indexing is necessary for pbaa
	
	tag "${sample}"
	// publishDir params.trimmed_fastq, pattern: "*.fai" pattern mode: 'copy'
	
	cpus 1
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	output:
	tuple path("*_trimmed.fastq"), val(sample), val(animal), emit: fastq
	path("*.fai"), emit: index
	
	script:
	"""
	mv ${fastq} ${sample}_trimmed.fastq
	samtools faidx ${sample}_trimmed.fastq.fai
	"""

}

process INDEX_GUIDE_FASTA {
	
	// make index from each potential guide FASTA file
	// indexing is necessary for pbaa
	
	punlishDir params.pbaa_resources, pattern: "*.fasta.fai", mode: 'copy'
	
	input:
	path(guide_fasta)
	
	output:
	path("*.fasta.fai")
	
	script:
	"""
	samtools faidx ${guide_fasta}
	"""

}

process RUN_PBAA {
	
	// run pbaa clustering on each FASTQ that has been trimmed
	// suppress stderr so you can see which samples fail
	
	tag "${sample}"
	publishDir params.pbaa_clusters, pattern: "*_passed_cluster_sequences.fasta", mode: 'copy'
	
	cpus 1
	
	input:
	tuple path(fastq), val(sample), val(animal)
	path(index)
	each path(guide_fastas)
	each path(guide_fasta_indices)
	
	output:
	path("*_passed_cluster_sequences.fasta"), val(sample), val(animal)
	
	script:
	if( sample.toLowerCase().contains("dpa") )
		"""
		pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		"MW679627.fasta" \
		${fastq} \
		"${sample}_DPA" 2> /dev/null
		"""
	else if( sample.toLowerCase().contains("dpb") )
		"""
		pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		"MW679628.fasta" \
		${fastq} \
		"${sample}_DPB" 2> /dev/null
		"""
	else if( sample.toLowerCase().contains("dqa") )
		"""
		pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		"MW679620.fasta" \
		${fastq} \
		"${sample}_DQA" 2> /dev/null
		"""
	else if( sample.toLowerCase().contains("dqb") )
		"""
		pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		"MW679622.fasta" \
		${fastq} \
		"${sample}_DQB" 2> /dev/null
		"""
	else if( sample.toLowerCase().contains("drb") )
		"""
		pbaa cluster \
		--min-read-qv 30 \
		--max-reads-per-guide 1000 \
		--max-alignments-per-read 2000 \
		--num-threads ${task.cpus} \
		"cy0333-drb.fasta" \
		${fastq} \
		"${sample}_DRB" 2> /dev/null
		"""

}

process RENAME_FOR_LAA {
	
	// pbaa adds extraneous sequence to the output fasta
	// rename to remove
	
	tag "${sample}"
	
	input:
	path(pbaa_cluster), val(sample), val(animal)
	
	output:
	path("*.fasta"), val(sample), val(animal)
	
	script:
	"""
	mv ${pbaa_cluster} ${sample}.fasta
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
	publishDir params.sample_clusters, mode: 'symlink'
	
	cpus 1
	
	input:
	path(fasta), val(sample), val(animal)
	
	output:
	path("*.fasta.gz"), val(sample), val(animal)
	
	script:
	"""
	gzip < /dev/null > ${sample}.fasta.gz && \
	dedupe.sh -Xmx1g ow=t \
	in=${sample}.fasta.gz \
	outbest=${sample}.fasta.gz \
	fo c \
	threads=${task.cpus}
	"""

}

process RENAME_CLUSTERS {
	
	// prepend sample name to cluster names. This makes it easier to track which sequences are associated with which clusters.
	
	tag "${sample}"
	
	input:
	path(fasta), val(sample), val(animal)
	
	output:
	path("*.fasta.gz"), val(sample), val(animal)
	
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
	
	tag "merging mamu"
	publishDir params.merged_clusters, mode: 'copy'
	
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
	
	tag "merging mafa"
	publishDir params.merged_clusters, mode: 'copy'
	
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
	
	input:
	tuple path(putative), val(animal)
	
	output:
	tuple path("*putative_alleles.fasta"), val(animal)
	
	script:
	"""
	#!/usr/bin/python
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	
	# parse input FASTA
	with open(animal + "_putative_alleles.fasta", "w") as handle:
		for idx, record in enumerate(SeqIO.parse(putative, "fasta")):
			record.id = str(current_time) + '-' + str(idx)
			SeqIO.write(record, handle, "fasta")
	"""

}

process PARSE_IPD_GENBANK {
	
	// use biopython to take IPD Genbank file and separate into FASTA
	// of genomic DNA (intron containing) sequences and cDNA sequences (no intron)
	
	tag "${animal}"
	publishDir params.ipd_ref_sep, mode: 'copy'
	
	cpus 1
	
	input:
	path(guide_fasta)
	
	output:
	tuple path("*gdna_reference.fasta"), val(animal_name), emit: gdna
	tuple path("*cdna_reference.fasta"), val(animal_name), emit: cdna
	
	script:
	animal_name = name.substring(8,12)
	
	"""
	ipd_to_gdna_cdna.py ${animal_name}
	"""

}

process MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA {
	
	// identify putative alleles whose sequences match full-length gDNA sequences already in IPD.
	// save BAM file of reads that map to known IPD alleles and FASTA.gz file of reads that do not.
	
	tag "${putative_animal}"
	
	cpus 1
	
	when:
	ref_animal == putative_animal
	
	input:
	each tuple path(gdna), val(ref_animal)
	tuple path(putative), val(putative_animal)
	
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
	
	when:
	ref_animal == putative_animal == fasta_animal
	
	input:
	tuple path(sam), val(putative_animal)
	tuple path(fasta), val(fasta_animal)
	each tuple path(gdna), val(ref_animal)
	
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
	
	when:
	putative_animal == ref_animal_name
	
	cpus workflow.cpus
	
	input:
	tuple path(cdna_matches), val(putative_animal)
	
	output:
	tuple path("*merged.aln"), val(putative_animal), emit: merged
	tuple path("*gdna_single_temp.fasta"), val(putative_animal)
	
	script:
	"""
	#!/usr/bin/python
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	
	# create gDNA sequence object
	for gdna_record in SeqIO.parse(input[0], "fasta"):
		# print status message
		print("Finding cDNA matches to " + str(gdna_record.name))
		
		# map gDNA to cDNA
		# returning FASTA file of cDNA within 50nt
		# write gDNA sequence to file
		with open(output[1], "w") as output_handle:
			SeqIO.write(gdna_record, output_handle, "fasta")
			
		# create cDNA sequence object from nearby matches
		for cdna_record in SeqIO.parse(input[1], "fasta"):
	
			# write a tab delimited line containing the gdna_record.id, cdna_record.id, cdna_record length, and match length
			# use command substitution to get the match length
			# run MUSCLE on sequence objects, but do it in a stream to avoid a lot of file I/O
			# pipe output to CLUSTALW format
			# then count the number of '*' characters that denote matches between the two sequences
			# this works for class I
			# for class II, the cDNA can be longer than the gDNA so this doesn't work
			# if the count of natching characters equals the number of cDNA characters, write to file
			# add maxiters = 2 to accelerate processing per https://www.drive5.com/muscle/manual/compromise.html
	
			shell( 'echo "' + gdna_record.name  + '\t' + cdna_record.name + '\t' + str(len(gdna_record.seq)) + '\t' + str(len(cdna_record.seq)) + '\t' + ' \
				$(echo ">' + gdna_record.name + '\n' + str(gdna_record.seq) + '\n>' + cdna_record.name + '\n' + str(cdna_record.seq) + '" \
				| muscle -maxiters 2 -quiet -clwstrict \
				| grep "^ " | grep "\*" -o | wc -l )" >> {output[0]}')
	"""
	

}

process FIND_MUSCLE_CDNA_GDNA_MATCHES {
	
	// run awk to find gdna sequences that fully match sequence of cdna
	
	tag "${putative_animal}"
	publishDir params.cdna_identical, mode: 'copy'
	
	input:
	tuple path(merged), val(putative_animal)
	
	output:
	tuple path("*matches.aln"), val(putative_animal)
	
	script:
	"""
	awk \'{{if( $4 == $5  ) print $0}}\' ${merged} > ${putative_animal}_matches.aln
	"""

}

process RENAME_MUSCLE_CDNA_MATCHES_FASTA {
	
	// input gDNA sequences
	// change the name of the sequence header
	// write output sequences that match cDNAs to new file
	
	tag "${gdna_animal}"
	publishDir params.cdna_identical, mode: 'copy'
	
	when:
	cdna_animal == gdna_animal
	
	input:
	tuple path(cdna_matches), val(cdna_animal)
	tuple path(matches_aln), val(gdna_animal)
	
	output:
	tuple path("*cdna_matches.fasta"), val(cdna_animal)
	
	script:
	"""
	#!/usr/bin/python
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	import csv
	
	# create output file in case it is empty
	shell('touch {output[0]}')
	
	# read input FASTA line-by-line
	for record in SeqIO.parse(input[0], "fasta"):
	
		# parse file with gdna sequences that match cdna sequences
		with open(input[1]) as tsvfile:
			reader = csv.reader(tsvfile, delimiter='\t')
			for row in reader:
	
				# test if name of sequence in cdna match file matches gdna sequence name
				if row[0] == record.name:
					# update name of sequence in output file
					record.description = row[1] + '|' + row[0]
	
					# write to file
					with open(output[0], "a") as handle:
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
	
	when:
	no-gdna_animal == cdna_animal
	
	input:
	tuple path(no-gdna_match), val(no-gdna_animal)
	tuple path(cdna_matches), val(cdna_animal)
	
	output:
	tuple path("*novel.fasta"), val(cdna_animal)
	
	script:
	if( cdna_matches.isEmpty() )
		"""
		touch ${no-gdna_animal}_novel.fasta
		"""
	else
		"""
		mapPacBio.sh in=${no-gdna_match} ref=${cdna_matches} outu=${no-gdna_animal}_novel.fasta subfilter=0
		"""

}

process MERGE_READS {
	
	// concatenate files to align with clustal omega
	
	tag "${novel_animal}"
	
	when:
	ref_animal == novel_animal == cdna_match_animal
	
	input:
	each tuple path(gdna_ref), val(ref_animal)
	tuple path(novel), val(novel_animal)
	tuple path(cdna_matches), val(cdna_match_animal)
	
	output:
	tuple path("*reads.fasta"), val(reads_animal)
	
	script:
	"""
	cat ${gdna_ref} ${novel} ${cdna_matches} > ${novel_animal}_reads.fasta
	"""

}

process CLUSTAL_ALIGN {
	
	// align reads with clustal omega
	// generate alignment in fasta format and distance matrix
	// distance matrix can be used to parse closest matches to known alleles
	
	tag "${novel_animal}"
	publishDir params.novel_alleles, pattern: '*distance.txt', mode: 'copy'
	
	when:
	reads_animal == novel_animal
	
	cpus workflow.cpus
	
	input:
	tuple path(reads), val(reads_animal)
	tuple path(novel), val(novel_animal)
	
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
	
	when:
	clustal_animal == novel_animal == match_animal
	
	input:
	tuple path(distances), val(clustal_animal)
	tuple path(novel), val(novel_animal)
	tuple path(cdna_matches), val(match_animal)
	
	output:
	tuple path("*novel_closest_matches.xlsx"), val(clustal_animal), emit: closest_matches
	tuple path("*distances_tmp.txt"), val(clustal_animal)
	
	script:
	"""
	parse_clustalo_distances.py ${clustal_animal}
	"""

}

process CREATE_GENOTYPING_FASTA {
	
	// to check accuracy of novel sequences, want to genotype IPD gDNA matches, cDNA matches, and novel sequences
	// need to create FASTA file that contains all of these putative_alleles along with their classification
	// 
	// 27319 - create FASTA file that only has cDNA extensions and novel alleles
	
	tag "${putative_animal}"
	publishDir params.genotyping, pattern: 'classified.fasta', mode: 'copy'
	publishDir params.putative_new, pattern: '*putative.fasta', mode: 'copy'
	
	when:
	putative_animal == ref_animal == match_animal == closest_animal
	
	input:
	tuple path(putative), val(putative_animal)
	each tuple path(gdna_ref), val(ref_animal)
	tuple path(matches), val(match_animal)
	tuple path(closest_matches), val(closest_animal)
	
	output:
	tuple path("*classified.fasta"), val(putative_animal), emit: classified
	tuple path("*putative.fasta"), val(putative_animal), emit: new_allele
	
	script:
	"""
	create_genotyping_fasta.py ${putative_animal}
	"""

}

process GENOTYPE_CSS {
	
	// identify reference sequences that are fully contained in CCS reads
	// use trimmed full-length reads for mapping
	// use higher minratio and smaller maxindel per emails with Brian Bushnell
	
	tag "${sample}"
	publishDir params.genotyping, pattern: '*.sam', mode: 'symlink'
	
	when:
	animal == classified_animal
	
	cpus = 1
	
	input:
	tuple path(fastq), val(sample), val(animal)
	each tuple path(classified), val(classified_animal)
	
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
	subfilter=0 \
	threads=1 \
	ow=t
	"""

}

process CREATE_MAMU_GENOTYPING_CSV {
	
	// create CSV file containing sample name, allele name, and count of reads matching allele
	
	input:
	path(sam_list)
	
	output:
	tuple path("*.csv"), val("mamu")
	
	script:
	"""
	create_genotyping_csv.py "mamu"
	"""

}


process CREATE_MAFA_GENOTYPING_CSV {
	
	// create CSV file containing sample name, allele name, and count of reads matching allele
	
	input:
	path(sam_list)
	
	output:
	tuple path("*.csv"), val("mafa")
	
	script:
	"""
	create_genotyping_csv.py "mafa"
	"""

}


process CREATE_GENOTYPING_PIVOT {
	
	// make Excel-formatted pivot table from genotyping CSV
	
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
	
	input:
	tuple path(fasta), val(animal)
	
	output:
	tuple path("*mapped.gff"), val(animal), emit: gff
	tuple path("*gdna_single_temp.fasta"), val(animal)
	
	script:
	"""
	#!/usr/bin/python
	
	from Bio import SeqIO
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	import csv
	
	# create output file in case it is empty
	if os.stat(fasta).st_size == 0:
		shell('touch fasta')
	else:
		# read input FASTA line-by-line
		for record in SeqIO.parse(input[0], "fasta"):
			if 'dpa' in (record.name).lower():
				mrna_reference = config['annotate']['dpa']['mrna_reference']
				cds_annotation = config['annotate']['dpa']['cds_annotation']
			elif 'dpb' in (record.name).lower():
				mrna_reference = config['annotate']['dpb']['mrna_reference']
				cds_annotation = config['annotate']['dpb']['cds_annotation']
			elif 'dqa' in (record.name).lower():
				mrna_reference = config['annotate']['dqa']['mrna_reference']
				cds_annotation = config['annotate']['dqa']['cds_annotation']
			elif 'dqb' in (record.name).lower():
				mrna_reference = config['annotate']['dqb']['mrna_reference']
				cds_annotation = config['annotate']['dqb']['cds_annotation']
			elif 'drb' in (record.name).lower():
				mrna_reference = config['annotate']['drb']['mrna_reference']
				cds_annotation = config['annotate']['drb']['cds_annotation']
			
			with open(output[1], "w") as output_handle:
				SeqIO.write(record, output_handle, "fasta")
			
			# run exonerate        
			shell('exonerate \
			--showtargetgff \
			--showalignment FALSE \
			--showvulgar FALSE \
			--model cdna2genome \
			--query ' + mrna_reference + ' \
			--target {output[1]} \
			--refine full \
			--annotation ' + cds_annotation + ' \
			>> {output[0]}')
	
	"""

}

process PRELIMINARY_EXONERATE_PROCESS_GFF_PUTATIVE {
	
	// prepare exonerate GFF output for Geneious
	// need to pass experiment and paths to needed tar.gz files
	// script is in Docker container
	
	tag "${animal}"
	
	input:
	tuple path(gff), val(animal)
	
	output:
	tuple path("*processed.gff"), val(animal)
	
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
	
	input:
	tuple path(gff), val(animal)
	
	output:
	path("*putative.gff")
	
	script:
	"""
	awk \'{{if ($3 ~ /cds/) print $1"\t"$2"\t""CDS""\t"$4,"\t"$5"\t"$6"\t"$7"\t"$8"\t""Name=CDS;ID=CDS" }}\' ${gff} >> ${animal}_putative.gff
	"""

}

