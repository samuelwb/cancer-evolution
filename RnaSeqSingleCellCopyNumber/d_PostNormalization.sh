#!/bin/bash

# After step c you can plot your single-cell copy number data in R and you will see some cells lacking
# copy number alterations--these are the normal cells. (We recommend also plotting bulk DNA-Seq-based copy number 
# data for corroboration.) You will also see some cells with copy number alterations matching your bulk DNA-Seq 
# copy number data--these are the cancer cells. Make a "cell_identity.txt" file (second argument to this script)
# patterned after the example file provided in this folder, indicating whether each single cell is cancer or
# normal. This script will then further normalize copy number based on these normal cells that are *within your
# experiment* and thus are better controls than the exogenously provided normal cells from step b (i.e. HMECs).

inFile=				# input file from previous step CancerAndNormal-Patient1_SingleCellCombined_TiroshCnv.txt (from c_RemoveNormal.sh)
identityFile=cell_identity.txt	# your cell_identity.txt file describing whether each of your cells is normal or cancer, as described above.
outFile=			# your final output file show single-cell copy number analysis for each cell

python PostNormalization.py $inFile $identityFile $outFile

