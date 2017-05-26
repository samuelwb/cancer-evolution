#!/bin/bash

# Author: Stephen Piccolo

# Run VarScan copyCaller on the VarScan copynumber output to do copy calling, re-centering around 0, etc.

# User-defind inputs
copyNumFile=		# output file from the previous step (i.e. "Sample1_varscan.copynumber") 
minimumCoverage=15	# minimum coverage for region consideration
outFile=		# output file name (i.e. "Sample1_varscan.copynumber.called")
outFileHomDel=		# output file name for homozygous deletions (i.e. "Sample1_varscan.copynumber.homdel")

java -Xmx8g -jar VarScan.v2.3.9.jar copyCaller $copyNumFile --min-coverage $minimumCoverage --output-file $outFile --output-homdel-file $outHomDel


