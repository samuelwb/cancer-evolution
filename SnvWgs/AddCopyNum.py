# Author: Samuel Brady

# Takes as input a VCF file of mutations and a folder containing
# copy number data for the same samples and outputs the VCF file
# with only copy number 2 region mutations included

import sys
import os

inMatrix = sys.argv[1]
copyNumFolder = sys.argv[2]
copyNumExtension = sys.argv[3]
copyNumPrefix = sys.argv[4]
outFile = sys.argv[5]

vafMatrixFile = open(inMatrix)
outputFile = open(outFile, "w")


# get list of desired CNV files
cnvFileList = [thisFile for thisFile in os.listdir(copyNumFolder) if thisFile.endswith(copyNumExtension) and thisFile.startswith(copyNumPrefix)]
cnvFileList.sort()

# make a dictionary of good CNV regions where chromDict -> sampleDict -> startChromRegion-endChromRegion
chromList = [str(chromNum) for chromNum in range(1, 23)]
chromList.extend(["X", "Y", "M"])
chromDict = {}

for thisChrom in chromList:
	chromDict[thisChrom] = {}

	for cnvFile in cnvFileList:
		chromDict[thisChrom][cnvFile] = []

# now go through each input CNV file and populate chromDict
for cnvFile in cnvFileList:
	print cnvFile
	copyNumFile = open(os.path.join(copyNumFolder, cnvFile))

	for line in copyNumFile:
		lineList = line.split("\t")
		lineList[-1] = lineList[-1].rstrip("\n")

		chromosome = lineList[0].replace("chr", "")
		start = lineList[1]
		end = lineList[2]
		absoluteCopy = lineList[4]

		# skip chromosome M information
		if chromosome == "MT" or chromosome == "M":
			continue

		# if this region is not different from normal (between 1.5 and 2.5 absolute copy) add to chromDict
		if float(absoluteCopy) >= 1.5 and float(absoluteCopy) <= 2.5:
			chromDict[chromosome][cnvFile].append(start + "-" + end)

	copyNumFile.close()

# go through input VAF matrix file and output only lines with mutations in copy-number 2 regions in EVERY sample
counter = 0

header = vafMatrixFile.readline().rstrip("\n").split("\t")
header.append("Log2Copy")

outputFile.write("\t".join(header) + "\n")

for line in vafMatrixFile:
	if counter % 1000 == 0:
		print("On line " + str(counter) + " of input file")
	counter += 1

	lineList = line.rstrip("\n").split("\t")

	# get chromosome and location of this mutation
	chromosome, position = lineList[0].split(":")

	chromosome = chromosome.replace("chr", "")

	if chromosome == "M" or chromosome == "MT":
		lineList.append("Neutral")

		outputFile.write("\t".join(lineList) + "\n")
		continue

	# for this mutation, see if it's in the chromDict (i.e. same copy as normal sample, absolute copy between 1.5 and 2.5) for each sample
	inGoodRegion = True

	for cnvFile in cnvFileList:	
		goodRegionList = chromDict[chromosome][cnvFile]

		inGoodRegionThisFile = False

		for goodRegion in goodRegionList:
			goodRegionStart, goodRegionEnd = goodRegion.split("-")
			
			if int(position) >= int(goodRegionStart) and int(position) <= int(goodRegionEnd): # we found that this mutation is in a good region
				inGoodRegionThisFile = True
				break # if we found that this is in a good region, we can be done for this file and move on to the next file

		if not inGoodRegionThisFile:
			inGoodRegion = False
			break
	
	# add whether this mutation is in a copy-neutral or altered region
	if inGoodRegion:
		lineList.append("Neutral")

		outputFile.write("\t".join(lineList) + "\n")
	else:
		lineList.append("Altered")

		outputFile.write("\t".join(lineList) + "\n")

vafMatrixFile.close()
outputFile.close()


