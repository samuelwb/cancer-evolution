#!/bin/bash

# Author: Samuel Brady

# This bash script adjusts VAFs from previous step (step f) based on normal
# contamination determined by copy number or other analysis.

inVafMatrix=							# input VAF txt file for somatic variants from step f
inMultiplierFile=SampleMultipliersBasedOnNormalContam.txt	# for each sample in your analysis, this txt file provides a multiplier which VAFs will by multiplied by to account for normal contamination; for example, if tumor purity is 80%, the multiplier will be 1.25; see example file of indicated name for format
outVafMatrix=							# output VAF txt file with VAFs multiplied

python AdjustVaf.py $inVafMatrix $inMultiplierFile $outVafMatrix

