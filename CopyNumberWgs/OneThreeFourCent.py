# Author: Samuel Brady

# Stretches a CNV segmented file so that peaks touch absolute copy 1, 3, and 4
# to correct for normal contamination-induced collapsing of CNV data towards 2

import sys

inFile = sys.argv[1]
makeValuesCenteredAroundPlusMinusWhat = float(sys.argv[2])
outFile = sys.argv[3]

# read inFile into a matrix, converting to absolute copy as you go
inMatrix = []

inputFile = open(inFile)
thisSampleBestCenteredList = [] # will contain a list of input CNV values

for line in inputFile:
	lineList = line.rstrip("\n").split("\t")

	log2copy = float(lineList[-1]) # get the log2 fold change from zero copy number, which is the last item in the list
	absoluteCopy = 2 ** (1 + log2copy)

	lineList[-1] = absoluteCopy
	inMatrix.append(lineList)
	thisSampleBestCenteredList.append(absoluteCopy)

inputFile.close()

###############################
###############################

# NOW THAT YOU HAVE A LIST CENTERED AROUND 2 (thisSampleBestCenteredList) TRY DIFFERENT
# MULTIPLIERS TO ENSURE OTHER VALUES CENTER AROUND 1, 3, AND 4 AS WELL
thisSampleBestAllAroundList = thisSampleBestCenteredList # thisSampleBestAllAround will start out with the best-centered list then improve

# Now compute how many values fall near 1, 3, and 4 in best-centered list
howManyNear1or3or4original = 0
pm = makeValuesCenteredAroundPlusMinusWhat # "pm" = plus-minus; shortened variable name

for item in thisSampleBestCenteredList: # see if each item is in range
	if item != "NA": # if item is not NA...
		if item >= 1-pm and item <= 1+pm: # check if it's near 1
			howManyNear1or3or4original += 1 # and if so, increment the counter
		elif item >= 3-pm and item <= 3+pm:
			howManyNear1or3or4original += 1
		elif item >= 4-pm and item <= 4+pm:
			howManyNear1or3or4original += 1

# Try multiply (spreading) samples out around the axis of 2 using 1.0 to 4.0 as multipliers
for numToShift in [float(j) / 100 for j in range(100, 400, 1)]:
	multipliedColumnValues = []

	# Make the multipleColumnValues list by multiplying each value in thisSampleBestCenteredList
	for index in range(0, len(thisSampleBestCenteredList)):
		if thisSampleBestCenteredList[index] != "NA": # if not NA, shift the value and append to multipliedColumnValues
			multipliedColumnValues.append((float(thisSampleBestCenteredList[index])-2)*numToShift + 2)
		else: # if it is NA, just append NA
			multipliedColumnValues.append(thisSampleBestCenteredList[index])

	# Now count how many values are around 1, 3, and 4 in the range specified by user in shifted values
	howManyNear1or3or4shifted = 0

	for item in multipliedColumnValues: # see if each item is in range
		if item != "NA": # if item is not NA...
			if item >= 1-pm and item <= 1+pm: # check if it's near 1
				howManyNear1or3or4shifted += 1 # and if so, increment the counter
			elif item >= 3-pm and item <= 3+pm:
				howManyNear1or3or4shifted += 1
			elif item >= 4-pm and item <= 4+pm:
				howManyNear1or3or4shifted += 1

	if (howManyNear1or3or4shifted > howManyNear1or3or4original):
		thisSampleBestAllAroundList = multipliedColumnValues
		print("\tFound a better-centered dataset around 1,3,4 - updating: new numAround134 is " + str(howManyNear1or3or4shifted))
		print("\t\thowManyNear1or3or4shifted = " + str(howManyNear1or3or4shifted))
		print("\t\thowManyNear1or3or4original = " + str(howManyNear1or3or4original))
		print("\t\tnumToShift = " + str(numToShift))
		howManyNear1or3or4original = howManyNear1or3or4shifted

# now thisSampleBestAllAroundList has the values we desire
print("The input file had " + str(len(thisSampleBestCenteredList)) + " values")
print("The output file will have " + str(len(thisSampleBestAllAroundList)) + " values")

# modify inMatrix so that it has the new absolute copy values 
counter = 0

for thisValue in thisSampleBestAllAroundList:
	inMatrix[counter][-1] = str(thisValue)
	counter += 1

# output modified inMatrix to file
outputFile = open(outFile, "w")

for line in inMatrix:
	outputFile.write("\t".join(line) + "\n")

outputFile.close()


