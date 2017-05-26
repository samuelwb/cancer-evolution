#!/bin/bash

# Author: Samuel Brady

# Remove zero CNA lines from pileup file since these confuse VarScan2

# User-defined inputs
inPileup=	# input pileup file name (i.e. "Sample1.pileup")
outFile=	# output pileup file name (i.e. "Sample1.nozero.pileup")
	
python -u RemoveVarscanZeroLines.py $inPileup $outFile


