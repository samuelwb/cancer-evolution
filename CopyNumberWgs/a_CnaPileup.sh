#!/bin/bash

# Author: Samuel Brady

# Make pileup files from aligned BAM file in order to call copy number alterations

# User-defined input file
referenceGenomeFile=		# input reference genome fasta file here (i.e. "/data/Genomes/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa")
inBam=				# input BWA-generated aligned BAM file here (i.e. "Sample1.bam")
outPileup=			# desired output file name (i.e. "Sample1.pileup")

samtools mpileup -f $referenceGenomeFile $inBam > $outPileup

