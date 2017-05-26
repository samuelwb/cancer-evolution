#!/bin/bash

# Author: Samuel Brady

# This bash script adds a column to the VAF matrix from previous steps (output 
# of step g) indicating whether each variant is copy-altered or copy-neutral.

inFile=			# input VAF matrix file (which was output from previous step; now has VAFs multiplied to adjust for normal contamination)
copyNumFolder=		# folder containing copy number data for these samples
copyNumExtension=	# extension of copy number (segmented) files for these samples (i.e. "varscan.copynumber.called.segmented"); see folder "CopyNumberWgs" in this repository
copyNumPrefix=		# prefix of copy number (segmented) files for these samples (i.e. "QuadroCent"); see folder "CopyNumberWgs" in this repository
outFile=		# output VAF matrix with copy number status (neutral or altered) indicated as an additional column

python -u AddCopyNum.py $inFile $copyNumFolder $copyNumExtension $copyNumPrefix $outFile

