# Nextflow pipeline for discovering new MHC Class II alleles in Rhesus and Cynomolgus Macaques

## Overview

These pipeline uses macaque MHC class II genomic DNA amplicon data to identify high-quality novel MHC alleles, which can then be submitted to the [Immuno Polymorphism Database](https://www.ebi.ac.uk/ipd/). Important updates in this version include:

- When sample name, class II locus, and macaque species are specified in the input file name, the workflow tailors most steps to each locus and animal
- It uses pbAA [(Pacific BioSciences Amplicon Analysis)](https://github.com/PacificBiosciences/pbAA) to cluster similar amplicons, which are likely to represent discrete class II alleles.
- It also requires that each putative novel allele is present in at least two individual animals. 

## Quick Start
If Docker and NextFlow are already installed on your system, and simply run this command in a working directory of your choice:

```
nextflow run dholab/MHC-II-allele-discovery \
-latest \
--sample_manifest /path/to/sample_manifest \
--bam_folder /path/to/folder/of/PacBio_CSS_bam_files 
```

This command will automatically pull the workflow from GitHub and run it. If you do not have Docker and NextFlow installed, or want to tweak any of the workflow's default configurations, proceed to the following sections.

## Detailed Instructions

NOTE: You will need Git installed on a POSIXct-compatible system (e.g. MacOS, Linux, [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install)) to proceed.

Start by running `git clone` to pull the workflow into a directory of your choice, like so:

```
git clone https://github.com/dholab/MHC-II-allele-discovery.git .
```

After the workflow bundle has downloaded, you may want to ensure that the workflow's scripts are executable with `chmod +x bin/*`, though this should already be done.

You will also need to install the Docker engine if you haven't already. The workflow pulls all software dependencies Docker Hub. As such, you will never permanently install that software on your system. To install Docker, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/).

### Nextflow Installation

This workflow was built with the [NextFlow](https://www.nextflow.io/) workflow manager (v22.0.4). We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 21.04.0.5552`, you are set and ready to proceed.

To run the workflow, simply change into the workflow directory and run the following in the shell terminal:

```
nextflow run main.nf
```

If the workflow crashes for any reason, simply resume it with:

```
nextflow run main.nf -resume
```

### Workflow Compute Resources

This workflow will run most efficiently when at least eight cores, 16 gigabytes of RAM, and 32 gigabytes of disk are available to it.

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

This brief narrative was written by [David H. O'Connor](https://dho.pathology.wisc.edu/).

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
