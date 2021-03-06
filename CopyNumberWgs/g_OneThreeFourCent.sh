#!/bin/bash

# Author: Samuel Brady

# Converts log2 fold change from 0 copy to absolute copy and multiplies
# CNA values around the 2-axis (diploid axis) to make copy-loss and copy-gain peaks
# appear at 1, 3, and 4 since due to normal contamination present in nearly all tumor
# samples the CNA peaks for these copy-loss and copy-gain regions tend to collapse towards
# 2 and not be right on 1, 3, or 4 as they should be.

inFile=		# input file generated by previous step (i.e. "ZeroCent_Sample1_varscan.copynumber.called.segmented")
plusMinus=0.1	# region in which CNA peaks will be measured around 1, 3, and 4 (i.e. 0.1 means 0.9-1.1, 2.9-3.1, and 3.9-4.1 will be analyzed)
outFile=	# output file (i.e. "QuadroCentAbsolute_Sample1_varscan.copynumber.called.segmented")

python OneThreeFourCent.py $inFile $plusMinus $outFile

