#!/bin/bash

# Author: Stephen Piccolo

# Run segmentation on CNA data from VarScan

# User-defined inputs
copyCalledFile=		# output from previous script (i.e. "Sample1_varscan.copynumber.called")
outFile=		# segmented output file (i.e. "Sample1_varscan.copynumber.called.segmented")

Rscript --vanilla SegmentCopyNumber.R $copyCalledFile $outFile

