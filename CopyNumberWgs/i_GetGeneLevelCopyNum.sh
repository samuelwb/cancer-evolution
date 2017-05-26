#!/bin/bash

# This optional script gets the copy number for each gene from segmented data
# and is generally used on the output of e_Segment.sh

inFolder=./		# folder containing segmented file
inFolderPrefix=		# prefix of files in inFolder to analyze
inFolderExtension=	# suffix of files in inFolder to analyze (i.e. "varscan.copynumber.called.segmented")
transcriptFile=		# transcript definition file (i.e. "refGene.txt") which can be downloaded from hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
outFile=		# output file name, which will be a text file matrix of genes in rows, samples in columns (one per file in inFolder matching inFolderPrefix and inFolderExtension)
			# and log2 fold change copy number for each gene for each sample

python -u GetGeneLevelCopyNum.py $inFolder $inFolderPrefix $inFolderExtension $transcriptFile $outFile
