
import sys

inFile = sys.argv[1]
identityFile = sys.argv[2]
outFile = sys.argv[3]

# get cell identities (cancer vs. normal) from identityFile; store in dict
idDict = {} # dict where cell -> cellType[cancer vs. normal]

idFile = open(identityFile)

for line in idFile:
	lineList = line.rstrip("\n").split("\t")

	cell = lineList[0].replace("_", "-").replace("featureCounts", "tpm")
	cellType = lineList[1]

	idDict[cell] = cellType # add to dict

idFile.close()

# read inFile into rnaMatrix
rnaMatrix = []

inputFile = open(inFile)

for line in inputFile:
	lineList = line.split("\t")
	lineList[-1] = lineList[-1].rstrip("\n")

	rnaMatrix.append(lineList)

inputFile.close()

# determine the names and number of unique cancer specimens there are by checking for the names after the "_" in the header
specimenSet = set()
header = rnaMatrix[0]

for headerIndex in range(1, len(header) - 1):
	sampleName = header[headerIndex]

	specimen = sampleName.split("_")[1]
	specimenSet.add(specimen)

specimenList = list(specimenSet)
specimenList.sort()

# go through each cancer specimen and normalize to its normal cells
for thisSpecimen in specimenList:
	# find out which columns belong to this sample; store indices in goodColumns
	goodColumnsNormal = []
	goodColumnsCancer = []
	goodColumnsAll = []

	for headerIndex in range(1, len(header) - 1):
		thisSample = header[headerIndex]
		thisSampleProper, thisSampleSpecimen = thisSample.split("_")

		if thisSampleSpecimen == thisSpecimen: # see if we are working on the right specimen first of all, and if so, figure out if it's cancer or normal
			goodColumnsAll.append(headerIndex)
			cancerOrNormal = idDict[thisSampleProper]

			if cancerOrNormal == "Cancer":
				goodColumnsCancer.append(headerIndex)
				rnaMatrix[0][headerIndex] = thisSpecimen + "Malignant^" + rnaMatrix[0][headerIndex] # rename column header to enable sorting later
			elif cancerOrNormal == "Normal":
				goodColumnsNormal.append(headerIndex)
				rnaMatrix[0][headerIndex] = thisSpecimen + "Healthy^" + rnaMatrix[0][headerIndex] 
			else:
				sys.exit("ERROR: Invalid cell type.")

	# get average for each gene window for each normal sample; store in dict where geneWindow -> normalAverage
	normalDict = {}

	for rowIndex in range(1, len(rnaMatrix)):
		thisRowGene = rnaMatrix[rowIndex][0]
		thisRowNormalValues = [rnaMatrix[rowIndex][goodColIndex] for goodColIndex in goodColumnsNormal] #thisRowNormalValues = rnaMatrix[rowIndex][1:(normalColumns + 1)]

		thisRowNormalValues = [float(a) for a in thisRowNormalValues]

		thisRowMean = sum(thisRowNormalValues) / len(thisRowNormalValues)

		normalDict[thisRowGene] = thisRowMean

	# go through each value in the matrix (for this specimen) and define according to average
	for rowIndex in range(1, len(rnaMatrix)):
		gene = rnaMatrix[rowIndex][0]

		normalMean = normalDict[gene]

		for colIndex in goodColumnsAll:
			rnaMatrix[rowIndex][colIndex] = float(rnaMatrix[rowIndex][colIndex]) - normalMean

# now sort rnaMatrix so that we have normal then cancer (pre-tx) followed by normal then cancer (post-tx)
rnaMatrix = zip(*rnaMatrix) # first transpose to make it easier to sort by sample
rnaMatrix.sort()

geneRow = rnaMatrix.pop()
rnaMatrix.insert(0, geneRow)

rnaMatrix = [list(a) for a in zip(*rnaMatrix)]

# remove prefix from sample names (prefix was used for sorting)
for headerIndex in range(1, len(rnaMatrix[0]) - 1):
	rnaMatrix[0][headerIndex] = rnaMatrix[0][headerIndex].split("^")[1]

# output to file
outputFile = open(outFile, "w")

for line in rnaMatrix:
	outputFile.write("\t".join([str(a) for a in line]) + "\n")

outputFile.close()

