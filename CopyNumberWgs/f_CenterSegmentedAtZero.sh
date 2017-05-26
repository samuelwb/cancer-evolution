#!/bin/bash

# Author: Samuel Brady

# Center segmented data at 0 (absolute copy 2) since it is usually a bit off from 0

inFolder=		# folder containing previous step output files (i.e. "./")
inExtension=		# extension of segmented files (i.e. "varscan.copynumber.called.segmented"); all files in inFolder with this extension will be zero-centered
plusMinus=0.1		# region around zero in which the script will maximize the number of segments (i.e. 0.1 will count the number of segments from -0.1 to 0.1
			# and will shift the values to maximize the number of segments in that region
outFile=ZeroCent_	# prefix of output files (the full output file name will be the input file name with this pre-pended

python -u CenterSegmentedAtZero.py $inFolder $inExtension $plusMinus $outFile

