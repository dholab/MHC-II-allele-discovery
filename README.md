# Nextflow pipeline for discovering new MHC Class II alleles in Rhesus and Cynomolgus Macaques

## Introduction

These pipeline uses macaque MHC class II genomic DNA amplicons to identify high-quality novel MHC alleles, which can then be submitted to the [Immuno Polymorphism Database](https://www.ebi.ac.uk/ipd/). Important updates in this version include:

- When sample name, class II locus, and macaque species are specified in the input file name, the workflow tailors most steps to each locus and animal
- It uses pbAA [(Pacific BioSciences Amplicon Analysis)](https://github.com/PacificBiosciences/pbAA) to cluster similar amplicons, which are likely to represent discrete class II alleles.
- It also requires that each putative novel allele is present in at least two individual animals. 

### Configuration

The `nextflow.config` file specifies:

- the path to a TSV file containing the names of demultiplexed CCS files in BAM format and sample names. The sample names need to include the locus (DPA, DPB, DQA, DQB, or DRB) in the sample name; this is how the workflow knows which parameters to apply during processing
- the path to the folder containing these BAM files
- full-length gDNA reference sequences in the sense orientation for each class II locus. For this analysis, I use class II sequences from CY0333
- parameters for primer identification and filtering. Each locus has a different size and uses different amplifciation primers. Specify these here.
- guide sequences for PacBio pbaa clustering. For DPA, DPB, DQA, and DQB, I use the same CY0333 sequence as for orienting reads. For DRB, I used a FASTA file with the DRB sequences from CY0333.
- the allele database used for classification. The current class II amplicons amplify exon 2-exon 4 and include flanking sequence. I don't trust the existing gDNA sequences in IPD. In Geneious, I concatenated the exon 2 and exon 4 sequences from the `ipd-mhc-mamu-2022-07-11` release and then removed any intron sequences. This means there will be no matches to gDNA sequences during classification and allows only for cDNA extensions and novel alleles. This will allow us to create a high quality database of exon 2-4 gDNA sequences from scratch
- exon 2-4 reference sequences from each class II locus from CY0333 and a text file specifying the length of these sequences. Used by exonerate to apply putative exon boundaries to putative alleles

I modified the older version of this workflow, written for `Snakemake`, to:

- Create FASTQ files from the BAM files received from the UW Biotechnology Center
- Apply primer trimming and minimum size filtering for multiple loci at the same time. This works by looking for the locus name in the sample name.
- Remove an optimization when looking for cDNA matches that reduces the search space for MUSCLE. This was causing false negative classifications. I also added output that makes it easier to observe progress.
- Modified the annotation step to use different annotation files depending on the locus of the putative allele

The workflow was then refactored in `NextFlow`, at which point the ability to handle Rhesus and Cynomolgus Macaques simultaneously was added.

Invoke the workflow with:

`nextflow run main.nf`

## Narrative description

Here is a brief narrative of the workflow for use in contract reports and manuscripts:

A demultiplexed BAM file containing CCS reads for each sample was provided by the UW Biotechnology Center. The BAM file was converted to FASTQ format using `samtools bam2fq` (v1.12). The FASTQ reads were oriented in the same polarity by comparison with a reference sequence using `vsearch --orient` (v2.17.0_linux_x86_64). The NCBI Genbank accession numbers for the reference sequence used for orienting FASTQ reads are: {DPA: MW679627, DPB: MW679628, DQA: MW679620, DQB: MW679622, DRB: MW679619}. Following orienting, PCR primer sequences were removed sequentially from the 5' and 3' ends of each read using `bbduk.sh` (v38.90). Unexpectedly short reads that could not generate sequences spanning exons 2-4 were also removed with `bbduk.sh`. The minimum length filters and primer sequences used for trimming were:

```
  dpa:
	minimum_length: '2300'
	forward_primers: 'GTRTGCTATGTATTTGCTTCATAGGG'
	reverse_primers: 'GGGTAAGAGGTTAGATATGGGAGT'
  dpb:
	minimum_length: '5700'
	forward_primers: 'GGATTAGGTGAGAGTGGYGCC'
	reverse_primers: 'TAGGTCATGTTTGTTCCCCGAC,CCTGAGTACTTGGGACCTCATG'
  dqa:
	minimum_length: '1900'
	forward_primers: 'GGTTCTTTCTTCCCCTGTTCTCC'
	reverse_primers: 'ATAGTTTCAGTCAGCCCTGGATG'
  dqb:
	minimum_length: '4000'
	forward_primers: 'AGARCGCCCTGATCCCTCTAA'
	reverse_primers: 'YACCAAAGTTAAGGCTTGGTTCT'
  drb:
	minimum_length: '3500'
	forward_primers: 'GCGTCGCTGTCAGTGTCTTC'
	reverse_primers: 'TTCCTCCTCCAGAAAAGCCTATGG,AACCATGCACTGATSATTTCTGGA'
```

After trimming, reads from each sample were clustered using the Pacific Biosciences `pbaa` tool (v1.0.1). For all of the loci except Mamu-DRB, the guide FASTA sequences for `pbaa` are the same as the sequences used as references for clustering. For -DRB, NCBI accessions MW679619, MW679618, MW679617, and MW679616 were all used as guide sequences. Clustering was performed with these parameters, which optimize for sensitive allele discovery: `--min-read-qv 30 --max-reads-per-guide 1000 --max-alignments-per-read 2000`.

Clusters that are shared between two or more independent samples are considered as putative MHC class II alleles. These putative alleles are classified using a cascading approach. All of the rhesus macaque MHC class II sequences in the NHP IPD database as of July 11, 2022 were downloaded. Exon 2-3-4 sequences were concatenated for every sequence in this database release. Few IPD class II sequences have introns and those that do could contain intronic artifacts, so restricting the analysis to exon sequences is sensible. These IPD sequences were exhaustively mapped to the putative alleles using `MUSCLE` (v3.8.1551) with the parameters `-maxiters 2 -clwstrict`. Putative alleles that fully contain exact IPD exon 2-3-4 sequences are labeled as genomic DNA extensions of existing cDNA alleles. Putative alleles that do not match existing IPD sequences are denoted as novel class II alleles. The nearest IPD matches to these novel alleles were determined with `clustalo` (v1.2.4) using the `--full` parameter.

To visualize the distribution of putative alleles between samples, the initial FASTQ files from each sample were mapped to the cDNA extension and novel alleles using minimap2 (v2.20-r1061) using the `map-hifi` preset. Mapped reads were then filtered to retain those that do not have any substitutions relative to the reference sequence using `reformat.sh` (BBMap version 38.90) using the `subfilter=0` parameter.

Preliminary exon annotations were applied to the putative alleles using `Exonerate` (v2.4.0) with the `cdna2genome` model and `--refine full` The annotations were based on exon 2-3-4 versions of the same reference sequences used for orienting reads. These annotations were manually validated by eye in Geneious Prime® 2022.2.1.

This workflow is run in a pipeline built for `nextflow` (v22.0.4).
