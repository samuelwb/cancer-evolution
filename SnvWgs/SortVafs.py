# Author: Samuel Brady

# Sorts VAFs within each change pattern by high to low VAF

import sys
from operator import itemgetter

inFile = sys.argv[1]
outFile = sys.argv[2]

# get matrix header
inputFile = open(inFile)

header = inputFile.readline().rstrip("\n").split("\t")

changePatternIndex = header.index("SamplesWithMutation")

# read the file into a matrix
inMatrix = []
inMatrix.append(header)

for line in inputFile:
	lineList = line.rstrip("\n").split("\t")
	
	inMatrix.append(lineList)

inputFile.close()

# get the set of all change patterns
changePatternSet = set()

for rowIndex in range(1, len(inMatrix)):
	changePattern = inMatrix[rowIndex][changePatternIndex]

	changePatternSet.add(changePattern)

changePatternList = list(changePatternSet) # change set to list for convenience
changePatternList.sort()

# go through each change pattern and sort the VAFs; add this to sortedMatrix
sortedMatrix = []
sortedMatrix.append(header)

for changePatt in changePatternList:
	print("Working on " + changePatt)

	# get a matrix with only rows that have this changePattern
	subMatrix = []

	for rowIndex in range(1, len(inMatrix)):
		thisRowChangePattern = inMatrix[rowIndex][changePatternIndex]

		if thisRowChangePattern == changePatt:
			subMatrix.append(inMatrix[rowIndex])

	# now sort the VAFs by each sample belonging to this change pattern
	refIndex = header.index("Ref")

	samples = changePatt.replace("Samples:", "").split("-") # get list of samples with the mutation in this change pattern

	if changePatt == "Samples:None" or changePatt == "Samples:NoCoverage": # if this change pattern is a bad class, just output the mutations unsorted to the sortedMatrix
		sortedMatrix.extend(subMatrix)
		
		continue

	if changePatt == "Samples:All":
		samples = [index - 1 for index in range(2, refIndex)]

	columns = tuple([int(a) + 1 for a in samples]) # convert sample IDs to integers than add 1 to get the column for that sample; this gives us a list of columns that have samples with this mutation

	subMatrix = sorted(subMatrix, key = itemgetter(columns[0])) #lambda x: columns)
	
	sortedMatrix.extend([list(a) for a in subMatrix])

# output to file
outputFile = open(outFile, "w")

for line in sortedMatrix:
	outputFile.write("\t".join(line) + "\n")
	
outputFile.close()


