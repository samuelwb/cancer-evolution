# Author: Samuel Brady

# Gets the pattern of samples in which a mutation is present; this
# assumes that the first sample is germline, so only subsequent samples are analyzed 

import sys
from operator import itemgetter

inFile = sys.argv[1]
vafThresh = float(sys.argv[2]) # minimum VAF to be considered present, i.e. 0.05
outFile = sys.argv[3]

# go through input file, get change pattern, and add change pattern to outMatrix, which has all previous info plus one new column with change pattern
inputFile = open(inFile)
outMatrix = []

header = inputFile.readline().rstrip("\n").split("\t")
header.append("SamplesWithMutation") # add new column with change pattern

outMatrix.append(header)

refIndex = header.index("Ref") # this minus 1 is the last sample column, which helps us know where the sample columns end

for line in inputFile:
	lineList = line.rstrip("\n").split("\t")

	samplesWithMutation = "Samples:"

	inAllSamples = True
	hasNA = False

	# go through each sample in order and add information about what samples have this mutation
	for colIndex in range(2, refIndex):
		thisVaf = lineList[colIndex]

		if thisVaf == "NA": # if the VAF is NA, label the whole row has having an NA
			hasNA = True 
			break
		else:
			thisVaf = float(thisVaf)

		if thisVaf >= vafThresh:
			samplesWithMutation += str(colIndex - 1) + "-" # make a string of the form "2-4-5" where each number indicates a sample with the mutation
		else:
			inAllSamples = False
		
	if inAllSamples:
		samplesWithMutation = "Samples:All" # label mutations present in all samples as present in all

	if samplesWithMutation.endswith("-"):
		samplesWithMutation = samplesWithMutation[:-1] # get rid of trailing hyphen

	if hasNA:
		samplesWithMutation = "Samples:NoCoverage"

	if samplesWithMutation == "Samples:":
		samplesWithMutation = "Samples:None"

	lineList.append(samplesWithMutation)

	outMatrix.append(lineList)

inputFile.close()

# sort by last column ("SamplesWithMutation")
outMatrix = sorted(outMatrix, key = itemgetter(header.index("SamplesWithMutation"), header.index("Effect"), header.index("SangerGene?")))

header = outMatrix.pop()
outMatrix.insert(0, header)

# output to file
outputFile = open(outFile, "w")

for line in outMatrix:
	outputFile.write("\t".join(line) + "\n")

outputFile.close()

