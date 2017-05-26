#!/bin/bash

# Author: Stephen Piccolo

# Run VarScan copynumber on pileup to ascertain copy number

# User-defined inputs
dataRatio=		# the approximate ratio of normal / tumor coverage (i.e. 0.5 if you had 30X normal and 60X tumor coverage)
normalPileup=		# pileup from germline normal WGS sample (i.e. "Germline.nozero.pileup")
tumorPileup=		# pileup from tumor WGS sample (i.e. "Sample1.nozero.pileup")
outFile=		# output file prefix (i.e. "Sample1_varscan"; VarScan will add the ".copynumber" suffix to this)

java -Xmx8g -jar VarScan.v2.3.9.jar copynumber $normalPileup $tumorPileup $outFile --data-ratio $dataRatio

