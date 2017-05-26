# Author: Samuel Brady

# Multiplies VAFs in 1st argument file by multipliers in 2nd argument file;
# output file is 3rd argument

import sys

inVafMatrix = sys.argv[1]
inMultiplierFile = sys.argv[2]
outFile = sys.argv[3]

# get multipliers for each sample; store in dict where sample -> multiplier
multiplierFile = open(inMultiplierFile)

multiplierDict = {}

for line in multiplierFile:
	lineList = line.rstrip("\n").split("\t")

	sample = lineList[0]
	multiplier = lineList[1]

	multiplierDict[sample] = float(multiplier)

multiplierFile.close()

# go through input file; multiply as appropriate by the multipliers and output to outFile
inputFile = open(inVafMatrix)
outputFile = open(outFile, "w")

header = inputFile.readline().rstrip("\n").split("\t") # get header

outputFile.write("\t".join(header) + "\n") # write header to output file; this part won't change

for line in inputFile:
	lineList = line.rstrip("\n").split("\t")

	outRow = []

	for colIndex in range(0, len(lineList)):
		headerThisColumn = header[colIndex]

		if headerThisColumn not in multiplierDict: # if this isn't a sample column, just output it, without multiplication, to the output row
			outRow.append(lineList[colIndex])
		else: # but if it is a sample column, multiply by multiplier
			vaf = lineList[colIndex]

			if vaf == "NA":
				outRow.append("NA")
			else:
				vaf = float(vaf)

				vafMultiplied = vaf * multiplierDict[headerThisColumn]

				outRow.append(str(vafMultiplied))

	outputFile.write("\t".join(outRow) + "\n")

inputFile.close()
outputFile.close()

