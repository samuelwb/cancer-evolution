#!/bin/bash

# Author: Samuel Brady

# This bash script adds normal cell single-cell expression (from human mammary epithelial cells) as a diploid control from copy number normalization in single-cell RNA-Seq data.
# There are 85 individual HMECs' expression in the included file. We recommend adding 4x the number of HMECs as single cells in your experiment. Thus if you have approximately 80 single cells,
# add 4 batches of HMECs for normalization.

cancerTpm=			# input file of single-cell RNA-Seq TPM data, with genes in rows and single cell samples in columns (tab-delimited text file)
hmecTpm=HMEC_Combined.tpm	# an input file provide by Brady et. al of single-cell RNA-Seq data from 85 individual human mammary epithelial cells, used as a normal copy (diploid) control generated using the Fluidigm C1 system; found in this folder; may be less suitable for non-cancer types, thus users are encouraged to seek out "normal sample" gene expression data for their cancer type if this file does not work well for their purpose
outFile1="HMEC_"$cancerTpm	# first output file with 85 (1 batch) of HMECs added
outFile2="HMEC2_"$cancerTpm	# second output file with 85*2 (2 batches) of HMECs added
outFile3="HMEC3_"$cancerTpm	# third output file with 85*3 (3 batches) of HMECs added
outFile4="HMEC4_"$cancerTpm	# fourth output file with 85*4 (4 batches) of HMECs added; this will be your true output file, while the others can be deleted

python MergeMatrixOnRowNames.py $hmecTpm $cancerTpm $outFile1
python MergeMatrixOnRowNames.py $hmecTpm $outFile1 $outFile2
python MergeMatrixOnRowNames.py $hmecTpm $outFile2 $outFile3
python MergeMatrixOnRowNames.py $hmecTpm $outFile3 $outFile4 # true output

