#!/bin/bash

# Author: Samuel Brady

inFile="HMEC4_[YourInputFileFromPreviousStep]"	# input file from previous step (step a) with your single-cell RNA-Seq data combined with HMEC or other normal single-cell RNA-Seq data
refGeneFile=					# refGene.txt from UCSC (Google "refGene.txt" and you will find this file); this tells gene order along chromosomes for averaging gene windows
threshold=1700 					# remove samples (single cells) with fewer than this number of genes
geneThresh=30 					# remove genes that express an aggregate log2(tpm + 1) value less than this
excludeY=1 					# for now, use 1 to exclude Y
outFileNonWindow=				# output txt file name for normalized expression values without windows; not the primary output
outFileCellCycleScores=				# output txt file name for cell cycle scores for each cell; not the primary output
outFile=					# output txt file name with single-cell inferred copy based on 101-window gene averages; primary output

python -u SingleCellCna.py $inFile $refGeneFile $threshold $geneThresh $excludeY $outFileNonWindow $outFileCellCycleScores $outFile

