# Author: Samuel Brady

# Gets a CNA matrix for a matrix of single-cell RNA-Seq data (tpm)
# using the basic methods from Tirosh et al. (melanoma paper, Science 2016)

import sys
import math
from operator import itemgetter

inFile = sys.argv[1]
refGeneFile = sys.argv[2]
threshold = sys.argv[3] # threshold for number of genes a sample must express to include that sample
geneThresh = sys.argv[4] # threshold for sum of (log2(tpm + 1)) that a gene must express to be included
excludeY = sys.argv[5] # 1 if yes, 0 if no
outFileNonWindow = sys.argv[6]
outFileCellCycle = sys.argv[7]
outFile = sys.argv[8]

# **** make a matrix where you have appropriately log'ed the tpm values
logMatrix = []

inputFile = open(inFile)
header = inputFile.readline().split("\t")
header[-1] = header[-1].rstrip("\n")
logMatrix.append(header)

for line in inputFile:
	lineList = line.split("\t")
	lineList[-1] = lineList[-1].rstrip("\n")

	for index in range(1, len(lineList)):
		tpmVal = float(lineList[index])

		logTpmValPlus1 = math.log((tpmVal/10.0) + 1.0, 2)

		lineList[index] = logTpmValPlus1

	logMatrix.append(lineList)

inputFile.close()

# **** remove samples with fewer than $threshold number of genes
# go through each sample and calculate number of expressed genes
expressedDict = {}  # dict where sampleName -> number of expressed genes

for colIndex in range(1, len(logMatrix[0])):
	numExpressed = 0

	for rowIndex in range(1, len(logMatrix)):
		expressionVal = float(logMatrix[rowIndex][colIndex])

		if expressionVal != 0.0:
			numExpressed += 1

	sampleID = logMatrix[0][colIndex]

	expressedDict[sampleID] = numExpressed

# convert expressedDict to a list where first column is sample and second sample is number of genes so we can sort it
expressedList = []

for thisSample in expressedDict:
	numGenes = expressedDict[thisSample]

	expressedList.append([thisSample, numGenes])

expressedList = sorted(expressedList, key=itemgetter(-1))

# make a histogram for the user that indicates number of genes in each bin (500-gene bins up to 12,000)
bins = [i for i in range(0, 12000, 500)]
binDict = {} # dict containing number of samples in each bin (0 indicates 0 to 499, 500 indicates 500 to 999, etc.)

for thisBin in bins:
	binDict[thisBin] = 0

# go through each sample and assign it to the bins
for thisSample in expressedDict:
	numGenes = expressedDict[thisSample]

	currentBin = 0
	for thisBin in bins:
		if numGenes < thisBin:
			break
		else:
			currentBin = thisBin

	# add to binDict
	binDict[currentBin] += 1
	
# convert binDict to list where [bin, numSamples] and sort by bin
binList = []

for thisBin in binDict:
	numSamples = binDict[thisBin]

	binList.append([thisBin, numSamples])

binList.sort()

# now show a histogram showing the number of expressed genes
print("\n*** HISTOGRAM FOR NUMBER OF EXPRESSED GENES ***\n")

for rowIndex in range(0, len(binList)):
	thisBin = binList[rowIndex][0]
	numSamples = binList[rowIndex][1]

	print("Bin " + str(thisBin) + ":\t" + "o"*numSamples)

# remove samples with fewer than $threshold number of genes
threshold = int(threshold)
cleanedMatrix = []
goodColumns = [] # list of columns to include (where samples have more than threshold num of genes)
goodColumns.append(0) # we need to keep the first column which has the gene symbols

for colIndex in range(1, len(logMatrix[0])):
	sampleName = logMatrix[0][colIndex]

	numGenes = expressedDict[sampleName]

	if numGenes >= threshold:
		goodColumns.append(colIndex)

for rowIndex in range(0, len(logMatrix)): # go through each row in logMatrix and add to cleanedMatrix only those columns that are in the goodColumns list
	keepList = [logMatrix[rowIndex][j] for j in goodColumns]

	cleanedMatrix.append(keepList)

print("\nRemoved all samples with fewer than " + str(threshold) + " genes")

# **** keep only good genes
# get a list of genes to keep

inputFile = open(inFile)
inputFile.readline()

genesToKeep = []
genesToExclude = []

for line in inputFile:
	lineList = line.split("\t")
	lineList[-1] = lineList[-1].rstrip("\n")

	tpmLineSum = 0

	for index in range(1, len(lineList)):
		tpmVal = float(lineList[index])

		logTpmValPlus1 = math.log(tpmVal + 1.0, 2)
		tpmLineSum += logTpmValPlus1

	if tpmLineSum >= float(geneThresh):
		genesToKeep.append(lineList[0])
	else:
		genesToExclude.append(lineList[0])

print("\nWith an aggregateThresh of " + str(geneThresh) + " we will keep " + str(len(genesToKeep)) + " genes")
print("We will exclude " + str(len(genesToExclude)) + " genes")

inputFile.close()

# make a doubleCleanedMatrix (cleaned out both bad samples and bad genes)
doubleCleanedMatrix = []
doubleCleanedMatrix.append(cleanedMatrix[0]) # add header

for rowIndex in range(1, len(cleanedMatrix)):
	geneSymbol = cleanedMatrix[rowIndex][0]

	if geneSymbol in genesToKeep: # if it's a gene to keep, add it to the doubleCleanedMatrix
		doubleCleanedMatrix.append(cleanedMatrix[rowIndex])

# **** mean-center each gene
for rowIndex in range(1, len(doubleCleanedMatrix)):
	tpmList = [] # list of expression values for this row
	for index in range(1, len(doubleCleanedMatrix[rowIndex])):
		tpmList.append(float(doubleCleanedMatrix[rowIndex][index]))
		
	# get average expression value
	tpmAverage = sum(tpmList) / len(tpmList)

	# subtract average from each tpm value
	for index in range(1, len(doubleCleanedMatrix[rowIndex])):
		doubleCleanedMatrix[rowIndex][index] = doubleCleanedMatrix[rowIndex][index] - tpmAverage

# **** sort genes by chromosome
# make a dictionary where gene -> chrom_startPosition
genePositionDict = {}

refGene = open(refGeneFile)

for line in refGene:
	lineList = line.split("\t")

	chrom = lineList[2].split("chr")[1]
	startPosition = lineList[4]
	gene = lineList[12]

	if gene not in genePositionDict: # if the gene's not in there, add to dict
		genePositionDict[gene] = chrom + "_" + startPosition

refGene.close()

# add two columns to doubleCleanedMatrix that contain the chromosome and startPosition of the gene
doubleCleanedMatrix[0].append(0) # give both chromosome (this one) and chromPosition (line below) a header value of 0 so they will be sorted at the top
doubleCleanedMatrix[0].append(0)

for rowIndex in range(1, len(doubleCleanedMatrix)):
	thisGene = doubleCleanedMatrix[rowIndex][0]

	if thisGene not in genePositionDict:
		doubleCleanedMatrix[rowIndex].append("BAD-CHROM")
		doubleCleanedMatrix[rowIndex].append("BAD-POSITION")
		continue

	genePositionList = genePositionDict[thisGene].split("_")

	if len(genePositionList) > 2: # this indicates weird chromosome names, so skip it
                doubleCleanedMatrix[rowIndex].append("BAD-CHROM")
                doubleCleanedMatrix[rowIndex].append("BAD-POSITION")
		continue

	thisGeneChrom = genePositionList[0]
	thisGeneStart = genePositionList[1]

	# exclude Y if desired
	if excludeY == "1" and thisGeneChrom == "Y":
                doubleCleanedMatrix[rowIndex].append("BAD-CHROM")
                doubleCleanedMatrix[rowIndex].append("BAD-POSITION")
		continue

	# add chromosome and start site for this gene to the row
	if thisGeneChrom == "X" or thisGeneChrom == "Y":
                doubleCleanedMatrix[rowIndex].append(thisGeneChrom)
                doubleCleanedMatrix[rowIndex].append(int(thisGeneStart))
	else:
		doubleCleanedMatrix[rowIndex].append(int(thisGeneChrom))
		doubleCleanedMatrix[rowIndex].append(int(thisGeneStart))

# sort by chromosome followed by startPosition
doubleCleanedMatrix = sorted(doubleCleanedMatrix, key = itemgetter(-2, -1))

# **** eliminate genes that have a "BAD-POSITION" for their chromosome since we cannot use genes whose location we don't know
chromSortedMatrix = []

numberOfBadPositionGenes = 0

for line in doubleCleanedMatrix:
	if "BAD-POSITION" not in line:
		chromSortedMatrix.append(line)
	else:
		numberOfBadPositionGenes += 1

print("\nThere were " + str(numberOfBadPositionGenes) + " genes in non-informative/unclear chromosome locations; we will omit these")

# **** mean-center by sample
for colIndex in range(1, len(chromSortedMatrix[0]) - 2):
	# get the mean for this column; first get a list of values
	sampleValues = []

	for rowIndex in range(1, len(chromSortedMatrix)):
		sampleValues.append(float(chromSortedMatrix[rowIndex][colIndex]))

	averageForSample = sum(sampleValues) / len(sampleValues)

	# subtract mean from each value
	for rowIndex in range(1, len(chromSortedMatrix)):
		chromSortedMatrix[rowIndex][colIndex] -= averageForSample

# output mean-centered data 
outFileWithoutWindow = open(outFileNonWindow, "w")

for line in chromSortedMatrix:
	outFileWithoutWindow.write("\t".join([str(a) for a in line[:-2]]) + "\n")


outFileWithoutWindow.close()

# output cell cycle scores for G1/S and G2/M; the genes expressed in these two phases were determined from Tirosh, Garraway, Science 2016 melanoma paper (supp. table)
G1S = ['RPA2', 'CLSPN', 'NASP', 'USP1', 'DTL', 'EXO1', 'RRM2', 'MSH2', 'MCM6', 'CDCA7', 'MCM2', 'SLBP', 'GMNN', 'CASP8AP2', 'RFC2', 'MCM4', 'CCNE2', 'DSCC1', 'ATAD2', 'HELLS', 'RRM1', 'E2F8', 'FEN1', 'POLD3', 'RAD51AP1', 'PRIM1', 'UNG', 'UBR7', 'RAD51', 'WDR76', 'TIPIN', 'BLM', 'GINS2', 'CDC6', 'BRIP1', 'TYMS', 'UHRF1', 'PCNA', 'CHAF1B', 'CDC45', 'POLA1']

G2M = ['CDCA8', 'CDC20', 'KIF2C', 'PSRC1', 'ANP32E', 'CKS1B', 'NUF2', 'NEK2', 'CENPF', 'LBR', 'CENPA', 'BUB1', 'CKAP2L', 'HJURP', 'SMC4', 'ECT2', 'TACC3', 'CENPE', 'CDC25C', 'HMMR', 'TTK', 'ANLN', 'CDCA2', 'CKS2', 'TUBB4B', 'CDK1', 'KIF20B', 'KIF11', 'MKI67', 'CKAP5', 'NCAPD2', 'CDCA3', 'CBX5', 'TMPO', 'GAS2L3', 'CKAP2', 'G2E3', 'DLGAP5', 'NUSAP1', 'CCNB2', 'KIF23', 'CTCF', 'FAM64A', 'AURKB', 'TOP2A', 'HN1', 'BIRC5', 'NDC80', 'TPX2', 'UBE2C', 'AURKA', 'MCM5', 'RANGAP1', 'GTSE1']

#	get rows for each of these genes
G1Srows = []
G2Mrows = []

for rowIndex in range(1, len(chromSortedMatrix)):
	thisGene = chromSortedMatrix[rowIndex][0]

	if thisGene in G1S:
		G1Srows.append(rowIndex)

	if thisGene in G2M:
		G2Mrows.append(rowIndex)

#	go through each sample and get scores for these two modules (G1S and G2M); output to file
outCycle = open(outFileCellCycle, "w")
outCycle.write("Sample\tG1S_Score\tG2M_Score\n")

for colIndex in range(1, len(chromSortedMatrix[0]) - 2):
	G1SforThisSample = [] # list of expression values for G1S genes
	G2MforThisSample = [] # list of expressoin values for G2M genes for this sample

	for rowIndex in G1Srows:
		G1SforThisSample.append(float(chromSortedMatrix[rowIndex][colIndex]))

	for rowIndex in G2Mrows:
		G2MforThisSample.append(float(chromSortedMatrix[rowIndex][colIndex]))

	G1Sscore = (sum(G1SforThisSample) / len(G1SforThisSample)) + 1
	G2Mscore = (sum(G2MforThisSample) / len(G2MforThisSample)) + 1

	G1Sscore = str(G1Sscore)
	G2Mscore = str(G2Mscore)

	outCycle.write(chromSortedMatrix[0][colIndex] + "\t" + G1Sscore + "\t" + G2Mscore + "\n")

outCycle.close()

# **** replace all values below -3 with -3 and all values above 3 with 3
for colIndex in range(1, len(chromSortedMatrix[0]) - 2):
	for rowIndex in range(1, len(chromSortedMatrix)):
		thisValue = chromSortedMatrix[rowIndex][colIndex]

		if float(thisValue) > 3.0:
			chromSortedMatrix[rowIndex][colIndex] = 3
		
		if float(thisValue) < -3.0:
			chromSortedMatrix[rowIndex][colIndex] = -3

# **** get a 101-window
# first find locations where each chromosome ends; put in dict where chrom -> endRow
chromEndDict = {}

for rowIndex in range(2, len(chromSortedMatrix)):
	currentChrom = chromSortedMatrix[rowIndex][-2]
	previousChrom = chromSortedMatrix[rowIndex - 1][-2]

	if currentChrom != previousChrom:
		chromEndDict[previousChrom] = rowIndex - 1

# convert chromEndDict to a range where chrom -> [start, end]
chromRangeDict = {}

chromRangeDict[1] = [1, chromEndDict[1]]

for chromIndex in range(2, 23):
	chromEnd = chromEndDict[chromIndex]
	chromStart = chromEndDict[chromIndex - 1] + 1

	chromRangeDict[chromIndex] = [chromStart, chromEnd]

chromRangeDict["X"] = [chromEndDict[22] + 1, len(chromSortedMatrix) - 1]

print("\nObtaining 101-window average values...")

# go through each value in chromSortedMatrix and get the 101-window (50 upstream, 50 downstream, and itself)
for colIndex in range(1, len(chromSortedMatrix[0]) - 2):
	print("\tWorking on sample " + str(colIndex))

	for rowIndex in range(1, len(chromSortedMatrix)):
		# figure out the range that needs to be averaged (row range)
		thisChromosome = chromSortedMatrix[rowIndex][-2]

		thisChromRange = chromRangeDict[thisChromosome]
		thisChromStart = thisChromRange[0]
		thisChromEnd = thisChromRange[1]

		rangeStart = ""
		rangeEnd = ""

		if rowIndex - thisChromStart < 50: # if we are near the beginning of the chromosome...
			rangeStart = thisChromStart
			rangeEnd = rowIndex + 50
		elif thisChromEnd - rowIndex < 50: # if we are near the end of the chromosome...
			rangeStart = rowIndex - 50
			rangeEnd = thisChromEnd
		else: # "normal" case - we are near the middle of the chromosome
			rangeStart = rowIndex - 50
			rangeEnd = rowIndex + 50

		# get a list of values in the rangeToAverage
		rangeValues = []

		for rangeIndex in range(rangeStart, rangeEnd + 1):
			rangeValues.append(chromSortedMatrix[rangeIndex][colIndex])

		# get average
		rangeAverage = sum(rangeValues) / len(rangeValues)

		# convert this value to the range value
		chromSortedMatrix[rowIndex][colIndex] = rangeAverage

# *** remove chromosome positions (last 2 columns)
for rowIndex in range(0, len(chromSortedMatrix)):
	chromSortedMatrix[rowIndex].pop()
#	chromSortedMatrix[rowIndex].pop()

chromSortedMatrix[0][-1] = "ChromosomeNumber"

# output to file
outputFile = open(outFile, "w")

for line in chromSortedMatrix:
	outputFile.write("\t".join([str(a) for a in line]) + "\n")

outputFile.close()

