#!/bin/bash

# Author: Samuel Brady

# This bash script gets the pattern of change for each mutation
# based on the output of the previous step (step h). It assumes
# that samples in the original VCF were ordered chronologically, with
# the germline sample first (i.e. germline, day0_cancer, day100_cancer, 
# day200_cancer, etc.). It adds this pattern of change as a new column
# indicating this for each variant. Mutation clusters are determined
# based on the pattern of change for subclone analysis.

inFile=		# input VAF matrix (which was output from previous step, step h)
vafThresh=0.05	# VAF threshold for considering a sample present (i.e. "0.05" would mean VAF 0.01 is "not present" while 0.12 is "present")
outFile=	# output VAF matrix with change pattern indicated

python ChangePattern.py $inFile $vafThresh $outFile

