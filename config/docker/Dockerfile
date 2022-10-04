# Use miniconda as a parent image
# dockerreg.chtc.wisc.edu/dhoconno/pacbio-pbaa:25837
FROM continuumio/miniconda3

# update conda
RUN conda update -n base -c defaults conda

# install mamba
RUN conda install mamba=0.13.0 -n base -c conda-forge

# install tools needed for gdna amplicon allele discovery workflow
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
RUN mamba install -y -c conda-forge -c bioconda nextflow=22.04
RUN mamba install -y -c conda-forge -c bioconda bbmap=38.90
RUN mamba install -y -c conda-forge -c bioconda vsearch=2.17.0
RUN mamba install -y -c conda-forge -c bioconda samtools=1.12
RUN mamba install -y -c conda-forge -c bioconda biopython=1.78
RUN mamba install -y -c conda-forge -c bioconda pandas=1.2.4
RUN mamba install -y -c conda-forge -c bioconda pysam=0.16.0.1
RUN mamba install -y -c conda-forge -c bioconda minimap2=2.20
RUN mamba install -y -c conda-forge -c bioconda pigz=2.3.4
RUN mamba install -y -c conda-forge -c bioconda unzip=6.0
RUN mamba install -y -c conda-forge -c bioconda picard=2.23.9
RUN mamba install -y -c conda-forge -c bioconda exonerate=2.4.0
RUN mamba install -y -c conda-forge -c bioconda muscle=3.8.1551
RUN mamba install -y -c conda-forge -c bioconda pyfasta=0.5.2
RUN mamba install -y -c conda-forge -c bioconda gffread=0.12.1
RUN mamba install -y -c conda-forge -c bioconda ipython=7.19.0
RUN mamba install -y -c conda-forge -c bioconda openpyxl=3.0.5
RUN mamba install -y -c conda-forge -c bioconda clustalo=1.2.4

# install tools needed for demultiplexing PacBio CCS and making clusters
# pbaa does not work with python 3
# install miniconda v2
RUN wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh
RUN bash ~/miniconda.sh -b -p /miniconda2

# install lima and pbaa
RUN /miniconda2/bin/conda install -y -c bioconda lima=2.0.0 
RUN /miniconda2/bin/conda install -y -c bioconda pbaa=1.0.1

# install reference sequences used for read mapping
# MHC class I
# aligned HLA-A, HLA-B, HLA-C, HLA-E with MUSCLE and saved open reading frames plus 20bp flanking sequence
COPY ref/HLA-A-B-C-E-20bp-pad.fasta /ref/HLA-A-B-C-E-20bp-pad.fasta

# add script to convert SAM file to use 'D' instead of 'N' for introns
# preserves compatibility with Geneious
COPY scripts/fix-cigar.awk /scripts/fix-cigar.awk

# install files necessary for exonerate processing to Geneious-compatible GFF
COPY scripts/21295-exonerate_gff_to_alignment_gff3.pl /scripts/21295-exonerate_gff_to_alignment_gff3.pl
COPY scripts/process-gff.sh /scripts/process-gff.sh

# install HLA cDNA and annotation files for exonerate
# the annotation file is used to define where the CDS is located
# ONE
COPY ref/HLA-A-mRNA.fasta /ref/HLA-A-mRNA.fasta
COPY ref/HLA-A-mRNA-annotation.txt /ref/HLA-A-mRNA-annotation.txt

# install FASTA guide file and index used by pbaa
COPY ref/24923-vsearch-07-mhc-I.fasta /ref/24923-vsearch-07-mhc-I.fasta
COPY ref/24923-vsearch-07-mhc-I.fasta.fai /ref/24923-vsearch-07-mhc-I.fasta.fai
