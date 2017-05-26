#!/bin/bash

# Author: Stephen Piccolo

# From RNA-Seq fastq reads (either paired- or single-end) align reads to the human genome
# and determine gene expression (TPM, RPKM, and feature counts)

set -o errexit

# User-defined inputs
f1=			# fastq file for read 1 (i.e. "Sample1_R1.fastq")
f2=			# fastq file for read 2 (i.e. "Sample1_R2.fastq"); for single-end reads enter NULL
outFilePrefix=		# prefix for output files (i.e. "Sample1_rsub", which will yield output files such as "Sample1_rsub.tpm", "Sample1_rsub.rpkm", etc.)
referenceGenomeFile=	# reference genome (i.e. "/data/Genomes/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa")
gtfFile=		# GTF file showing gene locations (i.e. "/data/Genomes/iGenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf")
dataDir=		# output data directory name (i.e. "RsubreadOutput")

# make output directory if not already made
mkdir -p $dataDir

# create complete output prefix (including output folder specified by user)
outFilePrefix=$dataDir/$outFilePrefix

# execute RSubread
Rscript --vanilla AlignRnaSeq.R $referenceGenomeFile $f1 $f2 $gtfFile $outFilePrefix



