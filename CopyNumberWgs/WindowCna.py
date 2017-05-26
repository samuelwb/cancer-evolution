# Author: Samuel Brady

# Takes the window average of CNV data based on the
# window size input as argument 2

import sys

inFile = sys.argv[1]
windowLength = int(sys.argv[2])
outFile = sys.argv[3]

inputFile = open(inFile)
outputFile = open(outFile, "w")

currentChrom = "1"
currentChromWindow = []

counter = 0

for line in inputFile:
	if counter % 100000 == 0:
		print("Counter is at " + str(counter))
	counter += 1

	lineList = line.rstrip("\n").split("\t")

	chrom = lineList[0]
	absoluteCopy = float(lineList[-1])

	# check if we're still on the same chromosome
	if chrom == currentChrom:
		currentChromWindow.append(absoluteCopy)

		if len(currentChromWindow) == windowLength: # check whether we have reached windowLength yet
			# get window average and add to the output file
                        windowAverage = sum(currentChromWindow) / float(len(currentChromWindow)) # get window average
                        lineList[-1] = str(windowAverage)

                        outputFile.write("\t".join(lineList) + "\n")

			# remove first item for next iteration
			currentChromWindow.pop(0) # remove first item for next iteration
	else:
		currentChrom = chrom
		currentChromWindow = []
		currentChromWindow.append(absoluteCopy)

inputFile.close()
outputFile.close()

