#!/bin/bash

# Author: Samuel Brady

# This bash script uses output from step d (d_MarkSomatic.sh) to get variant allele frequencies (VAFs)
# for each variant. VAFs for each sample are output in a combined txt file. Separate txt files
# are generated for germline and somatic variants. 

inFile=			# input VCF generated from step d (FreeBayes variant calls marked as somatic or not somatic)
sangerCensusFile=	# Sanger Cancer Gene Census (csv) file indicating cancer genes, which are sorted at the top of the output so you can see them easily (download from http://cancer.sanger.ac.uk/census/)
outFileSomatic=		# output txt file showing VAFs of somatic variants
outFileGerm=		# output txt file showing VAFs of germline variants

python -u MakeVafMatrixFromVcf.py $inFile $sangerCensusFile $outFileSomatic $outFileGerm

