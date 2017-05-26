#!/bin/bash

# Author: Samuel Brady

# Take segment window averages of CNA data

segmentedFile=		# file generated from previous step (i.e. "QuadroCentAbsolute_Sample1_varscan.copynumber.called.segmented")
windowLength=30		# number of segments on which to obtain window averages
outFile=		# output windowed file name (i.e. "Windowed_QuadroCentAbsolute_Sample1_varscan.copynumber.called.segmented")

python WindowCna.py $segmentedFile $windowLength $outFile

