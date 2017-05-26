# Author: Samuel Brady

# Removes pileup file rows with zeros as copy-number since these
# don't work with VarScan2 

import sys

inFile = sys.argv[1] # input pileup  file
outFile = sys.argv[2] # output pileup file without zeros 

inputFile = open(inFile)
outputFile = open(outFile, "w")

excludedCounter = 0

for line in inputFile:
	lineList = line.split("\t")
	lineList[-1] = lineList[-1].rstrip("\n")

	copyNumber = lineList[3]
	
	# if copy number is zero, increment counter; otherwise write line to output file
	if copyNumber == "0":
		excludedCounter += 1
	else:
		outputFile.write(line)

print("For file " + inFile + " there were " + str(excludedCounter) + " lines with zeros that were excluded")
	

inputFile.close()
outputFile.close()

