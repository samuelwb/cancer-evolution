#!/bin/bash

# Author: Samuel Brady

# This script looks for fusion transcripts based on input RNA-Seq fastq files
# using TopHat-Fusion. It represents step 1; step 2 is to run TopHat-Fusion Post.

# User-defined inputs
thisRnaFile_1=		# fastq paired-end file #1
thisRnaFile_2=		# fastq paired-end file #2
outFolder=		# name of output folder, which must begin with "tophat_" (i.e. "tophat_Sample1")
fusionDatabasePath=	# path to database files included in TopHat installation (i.e. "/data/Software/tophat-2.0.6.Linux_x86_64/TopHat_Fusion_Database_Files/hg19")

mkdir $outFolder

/path/to/tophat-2.0.6.Linux_x86_64/tophat -o $outFolder -p 8 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --max-intron-length 100000 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM $fusionDatabasePath $thisRnaFile_1 $thisRnaFile_2

