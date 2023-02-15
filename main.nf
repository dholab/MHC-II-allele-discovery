#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	// input channels
	ch_sample_manifest = Channel
		.fromPath( params.sample_manifest )
		.splitCsv( header: ['fastq', 'sample', 'animal'] )
		.map { row -> tuple( "${params.spades_contig_folder}/${row.fastq}", row.sample, row.animal ) }
	
	ch_ipd_ref = Channel
		.fromPath( params.ipd_refs )
	
	// Workflow steps
	
	MAP_FASTQ (
		ch_sample_manifest
	)

	TRIM_TO_PACBIO_AMPLICONS (
		MAP_FASTQ.out
	)

	FILTER_HARD_CLIPPED_AMPLICONS (
		TRIM_TO_PACBIO_AMPLICONS.out
	)
	
	// TRIM_FASTQ (
	// 	ORIENT_FASTQ.out
	// )

	// CLUSTER_PER_SAMPLE (
	// 	TRIM_FASTQ.out.fastq,
	// )
	
	// RENAME_CLUSTERS (
	// 	CLUSTER_PER_SAMPLE.out
	// )
	
	// MERGE_PER_MAMU_CLUSTERS (
	// 	RENAME_CLUSTERS.out
	// 		.filter { it[2] == "mamu" }
	// 		.map { fasta, sample, animal -> fasta }
	// 		.collect()
	// )
	
	// MERGE_PER_MAFA_CLUSTERS (
	// 	RENAME_CLUSTERS.out
	// 		.filter { it[2] == "mafa" }
	// 		.map { fasta, sample, animal -> fasta }
	// 		.collect()
	// )
	
	// SHARED_ANIMALS (
	// 	MERGE_PER_MAMU_CLUSTERS.out
	// 		.mix (
	// 			MERGE_PER_MAFA_CLUSTERS.out
	// 		)
	// )
	
	// RENAME_PUTATIVE_ALLELE_CLUSTERS (
	// 	SHARED_ANIMALS.out.putative
	// )
	
	// PARSE_IPD_GENBANK (
	// 	ch_ipd_ref
	// )
	
	// MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA (
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out,
	// 	PARSE_IPD_GENBANK.out.gdna
	// )
	
	// FILTER_EXACT_GDNA_MATCHES (
	// 	MAP_SHARED_CLUSTERS_TO_FULL_LENGTH_GDNA.out,
	// 	RENAME_PUTATIVE_ALLELE_CLUSTERS.out
	// 		.map { fasta, animal -> fasta },
	// 	PARSE_IPD_GENBANK.out.gdna
	// )
	
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
	
	// MERGE_READS (
	// 	PARSE_IPD_GENBANK.out.gdna,
	// 	EXTRACT_NOVEL_SEQUENCES.out,
	// 	RENAME_MUSCLE_CDNA_MATCHES_FASTA.out
	// 		.map { matches, animal -> matches }
	// )
	
	// CLUSTAL_ALIGN (
	// 	MERGE_READS.out,
	// 	EXTRACT_NOVEL_SEQUENCES.out
	// 		.map { novel, animal -> novel }
	// )
	
	// PARSE_DISTANCES (
	// 	CLUSTAL_ALIGN.out.distances,
	// 	EXTRACT_NOVEL_SEQUENCES.out
	// 		.map { novel, animal -> novel },
	// 	FILTER_EXACT_GDNA_MATCHES.out.cdna_matches
	// 		.map { matches, animal -> matches }
	// )
	
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

process MAP_FASTQ {

	// use minimap2 to map SPAdes contigs to combined DPA, DPB, DQA, DQB, DRB reference FASTA
	// this reference sequence contains 1kb of flanking sequence at end of each gene
	// this is needed to ensure that indels near the PacBio amplicon ends map correctly
	// after mapping, use samtools to sort the minimap2 output and convert to BAM file

	tag "${sample}"
	publishDir params.map_fastq, mode: 'copy'
	
	cpus 1
	memory '2.5 GB'
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fastq), val(sample), val(animal)
	
	output:
	tuple path("*.bam"), val(sample), val(animal)

	script:
	"""
	minimap2 -ax asm20 \
	${params.combined_reference} \
	${fastq} \
	--sam-hit-only --eqx \
	| samtools sort - \
	> ${sample}_mapped.bam
	"""
}

process TRIM_TO_PACBIO_AMPLICONS {
	
	// use samtools ampliconclip to trim mapped sequences to PacBio amplicon coordinates
	// this needs a BED file with the amplicon coordinates
	// the coordinates are relative to the 1kb extended reference sequences
	// when using ampliconclip need to take advantage of strand information
	// to trim the correct end of the sequence
	// 9100bp tolerance is needed because this is where PacBio amplicon is relative to 1kb extended DRB reference sequence
	
	tag "${sample}"
	publishDir params.trim_to_pacbio_amplicons, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(bam), val(sample), val(animal)
	
	output:
	tuple path("*.bam"), val(sample), val(animal)
	
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
	
	// for sequences to be alleles, they need to span the entire PacBio amplicon
	// after running ampliconclip, sequences that have hard clips at both ends are the ones we want
	// use a custom script to filter these sequences
	
	tag "${sample}"
	publishDir params.filter_hard_clipped_amplicons, mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(bam), val(sample), val(animal)
	
	output:
	tuple path("*.bam"), val(sample), val(animal)
	
	script:
	"""
	samtools index ${sample}_trimmed.bam
	
	python ${baseDir}/bin/filter_hard_clipped_ends.py ${sample}_trimmed.bam ${sample}_filtered.bam
	"""
	
}

process TRIM_FASTQ {
	
	// use bbduk to remove 5' and 3' primer sequences
	// failing to remove these sequences leads to artificial alleles with primer sequences at ends
	// need to trim to minlength or else short sequences choke pbaa
	// since sequences are oriented, only need to trim ends in one orientation
	// when using SPAdes contigs pre-process with bbduk to filter per-locus reads
	// concatenate all trimmed sequences into single file
	
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
	
	script:
	"""
	bbduk.sh -Xmx1g \
	in=${fastq} \
	ref=${params.dpa_amplicon} \
	outm=stdout.fastq \
	mcf=0.1 \
	| bbduk.sh -Xmx1g \
	int=f \
	in=stdin.fastq \
	literal=${params.dpa_forward_primers} \
	ktrim=l k=8 qin=33 editdistance=1 \
	minlength=${params.dpa_minimum_length} \
	out=stdout.fastq \
	| bbduk.sh -Xmx1g int=f in=stdin.fastq \
	literal=${params.dpa_reverse_primers} \
	ktrim=r k=8 qin=33 editdistance=1 \
	minlength=${params.dpa_minimum_length} \
	append=t \
	out=${sample}_trimmed.fastq

	bbduk.sh -Xmx1g \
	in=${fastq} \
	ref=${params.dpb_amplicon} \
	outm=stdout.fastq \
	mcf=0.1 \
	| bbduk.sh -Xmx1g \
	int=f \
	in=stdin.fastq \
	literal=${params.dpb_forward_primers} \
	ktrim=l k=8 qin=33 editdistance=1  \
	minlength=${params.dpb_minimum_length} \
	out=stdout.fastq \
	| bbduk.sh -Xmx1g int=f in=stdin.fastq \
	literal=${params.dpb_reverse_primers} \
	ktrim=r k=8 qin=33 editdistance=1 \
	minlength=${params.dpb_minimum_length} \
	append=t \
	out=${sample}_trimmed.fastq

	bbduk.sh -Xmx1g \
	in=${fastq} \
	ref=${params.dqa_amplicon} \
	outm=stdout.fastq \
	mcf=0.1 \
	| bbduk.sh -Xmx1g \
	int=f \
	in=stdin.fastq \
	literal=${params.dqa_forward_primers} \
	ktrim=l k=8 qin=33 editdistance=1  \
	minlength=${params.dqa_minimum_length} \
	out=stdout.fastq \
	| bbduk.sh -Xmx1g int=f in=stdin.fastq \
	literal=${params.dqa_reverse_primers} \
	ktrim=r k=8 qin=33 editdistance=1 \
	minlength=${params.dqa_minimum_length} \
	append=t \
	out=${sample}_trimmed.fastq

	bbduk.sh -Xmx1g \
	in=${fastq} \
	ref=${params.dqb_amplicon} \
	outm=stdout.fastq \
	mcf=0.1 \
	| bbduk.sh -Xmx1g \
	int=f \
	in=stdin.fastq \
	literal=${params.dqb_forward_primers} \
	ktrim=l k=8 qin=33 editdistance=1  \
	minlength=${params.dqb_minimum_length} \
	out=stdout.fastq \
	| bbduk.sh int=f in=stdin.fastq \
	literal=${params.dqb_reverse_primers} \
	ktrim=r k=8 qin=33 editdistance=1 \
	minlength=${params.dqb_minimum_length} \
	append=t \
	out=${sample}_trimmed.fastq
	
	bbduk.sh -Xmx1g \
	in=${fastq} \
	ref=${params.drb_amplicon} \
	outm=stdout.fastq \
	mcf=0.1 \
	| bbduk.sh -Xmx1g \
	int=f \
	in=stdin.fastq \
	literal=${params.drb_forward_primers} \
	ktrim=l k=8 qin=33 editdistance=1  \
	minlength=${params.drb_minimum_length} \
	out=stdout.fastq \
	| bbduk.sh int=f in=stdin.fastq \
	literal=${params.drb_reverse_primers} \
	ktrim=r k=8 qin=33 editdistance=1 \
	minlength=${params.drb_minimum_length} \
	append=t \
	out=${sample}_trimmed.fastq
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
	dedupe.sh -Xmx1g ow=t \
	in=${sample}_trimmed.fastq \
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


