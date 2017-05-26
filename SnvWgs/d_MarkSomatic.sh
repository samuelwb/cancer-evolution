#!/bin/bash

# Author: Brent Pedersen

# This bash script marks variants as somatic or not somatic using
# the "somatics" framework developed by Brent Pedersen. The "somatics"
# executable is also found in this folder. This assumes your
# germline sample is the first sample in the VCF.

inVcfFile=	# input VCF (quality-filtered) from previous step (step c)
outVcfFile=	# output VCF with variants marked as somatic or not somatic

./somatics $inVcfFile > $outVcfFile

