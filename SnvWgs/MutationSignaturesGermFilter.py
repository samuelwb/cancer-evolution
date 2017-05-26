# Author: Samuel Brady

# Gets mutation context for each mutation, given a matrix
# of mutations. It adds the mutation context to the matrix
# and summarizes the number of mutations with each signature
# in each change pattern.

import sys

inFile = sys.argv[1] # matrix of mutations
fastaFolder = sys.argv[2] # folder containing a fasta file for each chromosome (indicating the sequencing of each chromosome)
germVafFilter = float(sys.argv[3]) # the germline VAF for each mutation must be below this to be considered (i.e. 0.001)
outFile = sys.argv[4] # matrix of mutations with mutation context added as last column
outSummaryFile = sys.argv[5] # summary of mutation contexts for each change pattern in inFile

# ***********first define a function that gets the mutation context for a  mutation************
# This function returns the mutation context (1 base before, the base itself, and 1 base after)
# given an input chromosome and site on that chromosome; the folder containing chromosome
# fasta files (i.e. "17.fa", etc.) is also an input to the function

# Note that inputting the first or last base on a chromosome will break this function

import os

def getMutationContextFast(chromFastaFolder, chr, site):
        # open chromosome fasta file
        chrFileFast = open(os.path.join(chromFastaFolder, chr + ".fa"))
        chrFileFast.readline()

        # since we want to get the site before, let's start by getting site - 1
        siteBefore = int(site) - 1

        # determine the number of newlines that will be skipped to get to the site of interest
        numNewLinesToAdd = (siteBefore - 1) // 60

        # go to the site of interest
        chrFileFast.seek(siteBefore + numNewLinesToAdd - 1, 1)
        returnVal = chrFileFast.read(4) # read 4 characters because we want 3 characters total, and 1 of them may be a newline

        # remove any newline characters
        if "\n" in returnVal:
                returnVal = returnVal.replace("\n", "")
        else:
                returnVal = returnVal[:-1] # if no newline character get rid of the 4th character read

        chrFileFast.close()

        return returnVal

# read inFile into a matrix (inMatrix)
inMatrix = []

inputFile = open(inFile)
header = inputFile.readline().rstrip("\n").split("\t")
inMatrix.append(header)

for line in inputFile:
	lineList = line.rstrip("\n").split("\t")

	if lineList[1] == "NA":
		continue
	
	# test whether this germline VAF in this row is less than germVafFilter; if so, add to inMatrix
	germVaf = float(lineList[1])

	if germVaf < germVafFilter:
		inMatrix.append(lineList)

inputFile.close()

# go through each row of inMatrix (representing 1 mutation) and get the mutation context; add it as last column
inMatrix[0].append("MutContext") # add a column header indicating mutation context

header = inMatrix[0]

for rowIndex in range(1, len(inMatrix)):
	chrom, site = inMatrix[rowIndex][0].split(":") # get mutation chromosome and site

	# the next two commands clean up the chromosome name
	chrom = chrom.replace("chr", "")

	if chrom == "M":
		chrom == "MT"

	# get mutation ref and alt alleles
	ref = inMatrix[rowIndex][header.index("Ref")]
	alt = inMatrix[rowIndex][header.index("Alt")]

	# if this is an indel (i.e. either ref or alt is greater or less than 1 in length) do not attempt to get a mutation signature (indicate "NotApp-Indel")
	if len(ref) != 1 or len(alt) != 1:
		inMatrix[rowIndex].append("NotApp-Indel")
		continue

	# get trinucleotide context
	mutContext = getMutationContextFast(fastaFolder, chrom, site)

	# if the ref allele is G or A, get the complement of this as the mutation context; this narrows things down for us and gives us fewer mutation contexts (and fewer redundant ones)
	if ref == "G" or ref == "A": # in this case we need to convert everything
		complementDict = { "G": "C", "A": "T", "C": "G", "T": "A"}

		# convert the ref and alt alleles to their complement
		ref = complementDict[ref]
		alt = complementDict[alt]

		# convert the mutation context to its reverse complement
		mutContext = complementDict[mutContext[2]] + complementDict[mutContext[1]] + complementDict[mutContext[0]]

		# expand the mutation context to be not just the 3-letter "ACA" but something like "C>G_ACA", to indicate the full mutation context
		mutContext = ref + ">" + alt + "_" + mutContext
	else: # if we already are starting with C or T as ref, we don't need to do any complementing; we simply construct the final mutContext string with what we have
		mutContext = ref + ">" + alt + "_" + mutContext

	inMatrix[rowIndex].append(mutContext)

# ********* now summarize mutation signatures for each change pattern; start by making a list of possible mutation contexts (there are 96 total and range from "C>A_ACA" to "T>G_TTT")
possibleContexts = [] # will contain all possible contexts
sixPossibleChanges = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"] # for each of these 6 there will be 16 contexts, making 96 total
fourBases = ["A", "C", "G", "T"]

for thisChange in sixPossibleChanges: # this triple loop makes a list of all possible contexts - we start with the 6 possible changes
	for firstBase in fourBases: # then we loop through the possible first bases ("A", "T", "G", or "C")
		for lastBase in fourBases: # then we loop through the possible last bases
			mutationContextString = ""
			mutationContextString += thisChange + "_" # we start making the mutation context string by showing the change (i.e. "G>A_"
			mutationContextString += firstBase # then we add the first base - the base prior to the mutation
			mutationContextString += thisChange[0] # then we add the base itself, which we know to be the first character in the "thisChange" string
			mutationContextString += lastBase # then we add the last base - the base after the mutation

			possibleContexts.append(mutationContextString) # then we add this to the "possibleContexts" list

# great! now we have the list of possible mutation contexts - now we need to make a matrix where columns are change patterns, rows are mutation contexts, and
# the data are the number of mutations with a specific mutation context for each change pattern
# we'll start by making the empty matrix of mutation context data
mutationContextMatrix = []
changePatternSet = set()

for rowIndex in range(1, len(inMatrix)): # go through each row in change pattern out matrix and add the change pattern to the first row of mutationContextMatrix
	changePattern = inMatrix[rowIndex][header.index("SamplesWithMutation")]
	changePatternSet.add(changePattern)

# add each change pattern as a column header in first row
firstRow = []
firstRow.append("MutationContext") # give the first column header first though

for changePatt in changePatternSet:
	firstRow.append(changePatt)

mutationContextMatrix.append(firstRow) # give the mutationContextMatrix its first row with headers indicating change patterns

for thisContext in possibleContexts: # now go through each possible mutation context and make them the row headers - fill in the data with zeros also
	nextRow = []
	nextRow.append(thisContext) # row header is the mutation context

	nextRow.extend([int("0") for chgPattern in changePatternSet]) # for every change pattern put a zero into the matrix for this row

	mutationContextMatrix.append(nextRow) # now add this row to the mutationContextMatrix

# now that we have the matrix set up, we just need to fill it up with numbers;
# to do this, we'll loop through inMatrix (our general data table), and every
# time we find a valid mutation context (non-indel) we add data to the mutationContextMatrix
for rowIndex in range(1, len(inMatrix)):
	# get mutation context and change pattern for this row
	mutationContext = inMatrix[rowIndex][header.index("MutContext")]
	thisRowChangePattern = inMatrix[rowIndex][header.index("SamplesWithMutation")]

	if mutationContext != "NotApp-Indel": # if this is a valid mutation context (i.e. not a "C>CA" type indel or other unusable case)
		# if this is not an "NA", find the row and column index in the numerical mutationContextMatrix (that counts occurrences of contexts for each change pattern)
		# and increment it by 1
		rowIndexInMutationContextMatrix = possibleContexts.index(mutationContext) + 1
		colIndexInMutationContextMatrix = mutationContextMatrix[0].index(thisRowChangePattern)

		mutationContextMatrix[rowIndexInMutationContextMatrix][colIndexInMutationContextMatrix] += 1 # increment by 1

# write to output file
outputFile = open(outFile, "w")

for line in inMatrix:
	outputFile.write("\t".join(line) + "\n")

outputFile.close()

outSummaryFile = open(outSummaryFile, "w")

for line in mutationContextMatrix:
	outSummaryFile.write("\t".join([str(a) for a in line]) + "\n")

outSummaryFile.close()

