# Author: Samuel Brady

# Removes HMEC data from HMEC_Patient2 combined and orders samples (columns)
# according to input file

import sys

inRnaFile = sys.argv[1]
hmecPrefix = sys.argv[2]
outFile = sys.argv[3]

# read RNA input file into matrix
rnaMatrix = []

rnaFile = open(inRnaFile)

for line in rnaFile:
	lineList = line.split("\t")
	lineList[-1] = lineList[-1].rstrip("\n")

	rnaMatrix.append(lineList)

rnaFile.close()

# get list of good columns from header
header = rnaMatrix[0]
goodColumns = []
goodColumns.append(0)

for headerIndex in range(1, len(header)):
	thisSample = header[headerIndex]

	if not thisSample.startswith(hmecPrefix):
		goodColumns.append(headerIndex)

# output rnaMatrix to file with only goodColumns (those columns that are not HMECs)
outputFile = open(outFile, "w")

for line in rnaMatrix:
	outputFile.write("\t".join([str(line[index]) for index in goodColumns]) + "\n")

outputFile.close()

