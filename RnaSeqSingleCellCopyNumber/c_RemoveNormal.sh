#!/bin/bash

# Author: Samuel Brady

# This bash script removes the normal cells (HMECs in our case) spiked in at step a for copy number normalization

inFile=			# output file from previous step (step b); the file name that was the final argument to bash script b_SingleCellCna.sh
hmecPrefix=12342X	# prefix of normal cell sample names (column names) to remove; if you used the HMEC_Combined.tpm file we provided in this folder, use "12342X"
outFile=		# output file name of single-cell copy number WITHOUT normal cells (normal cells have been removed)

python RemoveNormal.py $inFile $hmecPrefix $outFile

