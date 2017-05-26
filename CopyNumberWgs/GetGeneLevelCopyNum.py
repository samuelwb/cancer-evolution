#This file generates a matrix (file) of gene-level CNA values
#from VarScan2 files (.varscan.copynumber.segmented extension). The output
#files from VarScan2 have CNA data in segmented format, rather than gene-level, format.

#Download "refGenes.txt" (hg19) from UCSC or other websites by Google-ing "refGenes.txt";
#this is an input file that helps determined what genes are found in chromosomal segments

#Output is a large matrix of rows being genes (from genes of interest file), columns being 
#samples, and values being CNV values.

import os
import sys
from utilities import naturalSort

# Four command-line arguments are needed, namely...
folderPath = sys.argv[1] # folder containing VarScan2 output files (CNV from exome data)
			# on which circular binary segmentation has been performed
inputFilePrefix = sys.argv[2]
inputFileExtension = sys.argv[3] # usually ".varscan.copynumber.segmented"; all files in folderPath with this extension will be analyzed
annotationFilePath = sys.argv[4] # annotation file (must be hg19) containing gene locations ("refGenes.txt")
outputPathAndFile = sys.argv[5] # FULL PATH of output matrix with patients' CNV data

# Move to directory specified by user and get list of directories in that folder
os.chdir(folderPath)

# Get list of files in the input directory that possess the extension of interest
fileList = [f for f in os.listdir(folderPath) if (os.path.isfile(folderPath + f)) and f.endswith(inputFileExtension) and f.startswith(inputFilePrefix)] 
fileList = naturalSort(fileList)

# Read in annotation file from annotationFilePath inputted by user;
# find locations of all genes in the file and store
# them in a dictionary where geneOfInterest -> chrX:11111-22222
chromList = [str(chromNum) for chromNum in range(1, 23)]
chromList.extend(["X", "Y", "MT"])

genesOfInterestLocations = {} # list of dictionaries where each dict is a chromosome from chromList and those dicts are geneSymbol:geneStart-geneEnd
for thisChrom in chromList:
	genesOfInterestLocations[thisChrom] = {}


uniqueGenes = set()

for line in file(annotationFilePath):
        # get file line and gene symbol for that line
	thisLine = line.split()
	geneSymbol = thisLine[12] # gene is in column 12
	uniqueGenes.add(geneSymbol)

        # if this line has a gene of interest, add its chromosomal location to dictionary
	chromosome = thisLine[2].split("chr")[1]
	geneStart = thisLine[4]
	geneEnd = thisLine[5]
		
	# if this is a weird chromosome then skip this line
	if chromosome not in genesOfInterestLocations:
		print("Skipping gene on chromosome " + chromosome)
		continue

	if geneSymbol not in genesOfInterestLocations:
		genesOfInterestLocations[chromosome][geneSymbol] = str(geneStart) + "-" + str(geneEnd)
	else: # if we have already seen this gene symbol, make the gene's coordinates shorter if this
                     # entry has shorter coordinates than previous entries
                     # (we want the smallest possible transcript length because some rare variants are extremely long)
		oldRange = genesOfInterestLocations[chromosome][geneSymbol]
		oldStart, oldEnd = oldRange.split("-")
		oldStart = int(oldStart)
		oldEnd = int(oldEnd)

		if int(geneStart) > oldStart:
			updatedStart = geneStart
		else:
			updatedStart = oldStart

		if int(geneEnd) < oldEnd:
			updatedEnd = geneEnd
		else:
			updatedEnd = oldEnd

		if updatedEnd > updatedStart: # assure start site is before end site
			genesOfInterestLocations[chromosome][geneSymbol] = str(updatedStart) + "-" + str(updatedEnd)

# outMatrix will hold the file output for later; for now we will give it a header of gene names found in refGene.txt
genesOfInterest = list(uniqueGenes)
genesOfInterest.sort()

#print("Length of dict is " + str(len(genesOfInterestLocations)))
#print("Length of set is " + str(len(uniqueGenes)))
#sys.exit("Made genesOfInterest dictionary!")

outMatrix = []
headerLine = []
headerLine.append(" ")
headerLine.extend(genesOfInterest)
outMatrix.append(headerLine)

# get CNV data for genesOfInterest for each patient file
for thisFileName in fileList:
	print("Converting " + thisFileName + " to gene-level CNV values")

	# go through file line-by-line; see if it has genesOfInterest
	# if so, add them to the thisSampleCNVdict (a dictionary where geneOfInterest -> CNV level)
	thisSampleCNVdict = {}
	thisSampleFile = file(thisFileName)

	counter = 0

	# look through each line of the file and see if those coordinates contain genes of interest
	for line in thisSampleFile:
		counter += 1
		
		if counter % 100 == 0:
			print("\tProcessing line " + str(counter) + " of " + thisFileName) 
		
		# retrieve chromosome segment information for this line
		thisLine = line.split()
		thisLine[-1] = thisLine[-1].rstrip("\n")
		
		chromosome = thisLine[0]
		startSite = int(thisLine[1])
		endSite = int(thisLine[2])
		segmentMean = thisLine[4]

		# check whether each geneOfInterest on this chromosome is in this chromosome segment
		geneDictThisChrom = genesOfInterestLocations[chromosome.replace("chr", "")]
		
		for gene in geneDictThisChrom:
			thisGeneSite = geneDictThisChrom[gene]
			thisGeneStart, thisGeneEnd = thisGeneSite.split("-") 
			thisGeneStart = int(thisGeneStart)
			thisGeneEnd = int(thisGeneEnd)

			#if this gene is in the range add it to the gene->segmentMean dictionary
			if thisGeneStart >= startSite and thisGeneStart <= endSite:
				thisSampleCNVdict[gene] = segmentMean

        # make a list of the genes of interest and their CNV levels
	thisSampleCNV = []
	thisSampleCNV.append(thisFileName) #first item is sample ID (fileName ID)

        # add CNV levels for each gene
	for thisGene in genesOfInterest:
		if thisGene in thisSampleCNVdict:
			thisSampleCNV.append(thisSampleCNVdict[thisGene])
		else:   #if gene does not appear in the CNVdict, just assign it the value of "NA"
			thisSampleCNV.append("NA")
			
	# add this patient's information to a combined cancer matrix for later output (outMatrix)
	outMatrix.append(thisSampleCNV)


# now transpose outMatrix to make patients in columns, genes in rows
outMatrix = zip(*outMatrix)

# output the outMatrix to a file
outTxtFile = file(outputPathAndFile, "w")

for line in outMatrix:
        outTxtFile.write("\t".join(line) + "\n")

