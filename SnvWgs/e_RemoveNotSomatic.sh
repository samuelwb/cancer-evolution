#!/bin/bash

# Author: Samuel Brady

# This bash script removes variants marked as not somatic (germline) by the previous step,
# step d.

inFile=		# VCF from previous step (with somatic status marked)
outFile=	# VCF with not somatic (germline) variants removed to get a somatic-only VCF

python RemoveNotSomatic.py $inFile $outFile

