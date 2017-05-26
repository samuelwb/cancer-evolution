#!/bin/bash

# Author: Samuel Brady

# This bash script removes variant scores with quality below 5 (or other
# value). Please modify the "QUAL" score below if a threshold other than
# 5 is desired.

inFile=		# VCF file of FreeBayes output annotated by SnpEff (from previous step, step b)
outFile=	# output VCF file with low-quality variants removed

/data/Software/freebayes/vcflib/bin/vcffilter -f "QUAL > 5" $inFile > $outFile 

