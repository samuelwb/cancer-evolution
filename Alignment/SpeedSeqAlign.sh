#!/bin/bash

# Author: Samuel Brady

# This script runs BWA alignment on input fastq files (i.e. from WGS) using the Speed-Seq suite

# User-defined inputs
referenceGenome=	# input your reference genome fasta file here (i.e. "/data/Genomes/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/BWAIndex/genome.fa")
fastqFile1=		# input fastq file pair 1 here (i.e. "Sample1_R1.fastq")
fastqFile2=		# input fastq file pair 2 here (i.e. "Sample1_R2.fastq")
prefix=			# input fastq prefix here for output file generation (i.e. "Sample1")

# Run BWA with 5 threads and 20 GB memory
speedseq align -o $prefix -M 20 -t 5 -R "@RG\tID:"$prefix"\tSM:"$prefix"\tLB:lib1" $referenceGenome $fastqFile1 $fastqFile2

