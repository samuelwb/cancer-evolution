# Author: Samuel Brady

# Removes germline mutations from VCF file

import sys

inFile = sys.argv[1]
outFile = sys.argv[2]


inputFile = open(inFile)
outputFile = open(outFile, "w")

for line in inputFile:
	if "NOT_SOMATIC" not in line:
		outputFile.write(line)

inputFile.close()
outputFile.close()
