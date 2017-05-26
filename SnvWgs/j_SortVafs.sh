#!/bin/bash

# Author: Samuel Brady

# This bash script sorts the VAF matrix output from the
# previous step (step i) based on change pattern for 
# subsequent attractive heatmap output.

inMatrix=	# input VAF matrix (from previous step, step i)
outMatrix=	# output VAF matrix (sorted)

python SortVafs.py $inMatrix $outMatrix

