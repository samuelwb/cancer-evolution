# Author: Samuel Brady

# Centers segmented CNV data around zero

import sys
import os

inFolder = sys.argv[1]
inExtension = sys.argv[2]
makeValuesCenteredAroundPlusMinusWhat = float(sys.argv[3])
outFilePrefix = sys.argv[4]

# get list of files
fileList = [thisFile for thisFile in os.listdir(inFolder) if thisFile.endswith(inExtension)] 

# go column-by-column (sample-by-sample) and normalize to near-integer values
for thisFile in fileList:
        print("\nNow centering file " + thisFile + " around 0")

	inFile = open(thisFile)

        # get the values for this file
        thisColumnValues = []

	counter = 0
        for line in inFile:
		if counter % 100000 == 0:
			print("\tReading in line " + str(counter))
		counter += 1

		lineList = line.split("\t")
		lineList[-1] = lineList[-1].rstrip("\n")
		cnvValue = lineList[-1]

		thisColumnValues.append(float(cnvValue))

	inFile.close()

        # now compute how many values fall within the range indicated by the user
        howManyNear2original = 0
        bottomOfRange = 0 - makeValuesCenteredAroundPlusMinusWhat # if input = 0.2, this will be 1.8
        topOfRange = 0 + makeValuesCenteredAroundPlusMinusWhat # if input = 0.2 this will be 2.2

	counter = 0
        for item in thisColumnValues: # see if each item is in range
                if counter % 100000 == 0:
			print("\t\tCounting line " + str(counter))
		counter += 1

		if item != "NA": # if item is not NA...
                        if item >= bottomOfRange and item <= topOfRange: # check if it's within the input range
                                howManyNear2original += 1 # and if so, increment the counter

        # this variable will hold the output data; every time a better output data list (better centered around 2)
        # is found this variable will be reassigned to that list
        thisSampleBestCenteredList = thisColumnValues # initialize it with thisColumnValues; then see if you can find something better

        # now center the values around 2 by adding or substracting increments of 0.01 (from -0.7 to 0.7) and see which centers the data around 2 best
        # (i.e. which gives the most values in the range input by the user surrounding 2.
	for numToShift in [float(j) / 100 for j in range(-70, 70, 1)]:
		print("\t\t\tOn numToShift " + str(numToShift))
                shiftedColumnValues = [] # This will hold this column's values shifted by numToShift

                # make the shiftedColumnValues list by shifting each value in thisColumnValues
                for index in range(0, len(thisColumnValues)):
                        if thisColumnValues[index] != "NA": # if not NA, shift the value and append to shiftedColumnValues
                                shiftedColumnValues.append(float(thisColumnValues[index]) + numToShift)
                        else: # if it is NA, just append NA
                                shiftedColumnValues.append(thisColumnValues[index])

                # now count how many values are around 2 in the range specified by user in shifted values
                howManyNear2shifted = 0

                for item in shiftedColumnValues: # see if each item is in range
                        if item != "NA": # if item is not NA...
                                if item >= bottomOfRange and item <= topOfRange: # check if it's within the input range
                                        howManyNear2shifted += 1 # and if so, increment the counter

                # if shifted values have more values centered around 2, make this shifted list the new "best list"
                if (howManyNear2shifted > howManyNear2original):
                        thisSampleBestCenteredList = shiftedColumnValues
                        print("\tFound a better-centered dataset around 2 - updating: new numAround2 is " + str(howManyNear2shifted))
                        print("\t\thowManyNear2shifted = " + str(howManyNear2shifted))
                        print("\t\thowManyNear2original = " + str(howManyNear2original))
                        howManyNear2original = howManyNear2shifted

	# now output to output file the same data as inFile except the last column (CNV information) changed to thisSampleBestCenteredList
	counter = 0

	outputFile = open(outFilePrefix + thisFile, "w")

	inFile = open(thisFile)

	for line in inFile:
		lineList = line.split("\t")
		lineList[-1] = lineList[-1].rstrip("\n")
		print counter

		lineList[-1] = thisSampleBestCenteredList[counter]

		outputFile.write("\t".join([str(a) for a in lineList]) + "\n")

		counter += 1
	
	outputFile.close()

	print("Finished on line " + str(counter) + " of file " + thisFile)
	print("The length of thisSampleBestCenteredList was " + str(len(thisSampleBestCenteredList)))
