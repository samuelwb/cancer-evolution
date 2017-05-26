# Author: Samuel Brady

# Converts a VCF file input to a matrix of variant allele
# frequencies (VAFs); also tells you whether each mutation touches
# a Sanger cancer gene census gene

import sys
from operator import itemgetter

inFile = sys.argv[1]
sangerCensusFile = sys.argv[2]
outFileSomatic = sys.argv[3]
outFileGerm = sys.argv[4]

# get list of sanger genes
sangerList = []

sangerFile = open(sangerCensusFile)
sangerFile.readline()

for line in sangerFile:
	sangerGene = line.split(",")[0]

	sangerList.append(sangerGene)

sangerFile.close()

# go through input VCF file; as you go, get VAF data and put into vafMatrixSomatic
inputFile = open(inFile)
vafMatrixSomatic = []
vafMatrixGerm = []

header = [] # this will hold the VCF header in list form
formatIndex = "" # this will hold the index of the "FORMAT" column in the VCF header

counter = 0

for line in inputFile:
	if counter % 10000 == 0:
		print("On line " + str(counter) + " of VCF")

	counter += 1

	if line.startswith("##"): # skip VCF header lines
		continue 
	elif line.startswith("#CHROM"): # get header information
		header = line.rstrip("\n").split("\t")

		# from this VCF file header, make a header appropriate for the output file
		outHeader = []
		
		outHeader.append("Site")

		formatIndex = header.index("FORMAT")	
		sampleList = header[(formatIndex + 1):] # get list of samples to add to outHeader

		outHeader.extend(sampleList)

		outHeader.append("Ref")
		outHeader.append("Alt")
		outHeader.append("Info")
		outHeader.append("Type")
		outHeader.append("Effect")
		outHeader.append("Gene")
		outHeader.append("ProteinChange")
		outHeader.append("SangerGene?")
		outHeader.append("Germline?")

		vafMatrixSomatic.append(outHeader)
		vafMatrixGerm.append(outHeader)
	else: # we have a legitimate row of data for a mutation; get the VAFs and output to outFileSomatic
		lineList = line.rstrip("\n").split("\t")
	
		outRow = [] # this will contain the output, including the VAFs

		# get chromosome site and position; add it as first entry to outRow
		chrom = lineList[0].replace("chr","")
		position = lineList[1]

		ref = lineList[3]
		alt = lineList[4]

		snpEffInfo = lineList[formatIndex - 1] # parse through snpEff info to get relevant information about the mutation (gene affected, missense, etc.)
		snpEffList = snpEffInfo.split("|")

		mutationType = snpEffList[1]
		mutationImpact = snpEffList[2]
		mutationGene = snpEffList[3]
		mutationProteinChange = snpEffList[10]

		if mutationType == "":
			mutationType = "NoMutationTypeListed"

		if mutationGene == "":
			mutationGene = "NoGeneListed"

		somaticNess = lineList[header.index("FILTER")]

		if mutationImpact == "MODERATE":
			mutationImpact = "INTERMEDIATE" # enables alphabetical sorting of variants by impact (HIGH, INTERMEDIATE, MODIFIER, LOW)

		if "," in alt: # if there are cases where there are multiple alternate loci in on VCF record, split it into multiple VCF records, as it were
			altList = alt.split(",")

			# get variant alelle frequencies - first need to get FORMAT information
			formatInfo = lineList[header.index("FORMAT")].split(":")

			aoIndex = formatInfo.index("AO")
			roIndex = formatInfo.index("RO")

			for altIndex in range(0,  len(altList)): # go through each alternate allele and get variant allele frequencies
				thisAlt = altList[altIndex]

				outRow = []

				outRow.append(chrom + ":" + position)

				for sampleIndex in range((formatIndex + 1), len(lineList)): # go through each sample for this alternate allele
					sampleRecord = lineList[sampleIndex].split(":")
					
					AOlist = sampleRecord[aoIndex].split(",")
					RO = sampleRecord[roIndex]

					VAF = ""

					if AOlist == "" or RO == "" or AOlist == "." or RO == ".":
						VAF = "NA"
					else:
						AO = AOlist[altIndex]

						if AO == "0" and RO == "0":
							VAF = "NA"
						else:
							VAF = float(AO) / (float(RO) + sum([float(a) for a in AOlist]))
				
					outRow.append(str(VAF))
		
				outRow.append(ref)
				outRow.append(thisAlt)
				outRow.append(snpEffInfo)
				outRow.append(mutationType)
				outRow.append(mutationImpact)
				outRow.append(mutationGene)
			
				if mutationProteinChange == "":
					outRow.append("None")
				else:
					outRow.append(mutationProteinChange)

				if mutationGene in sangerList:
					outRow.append("Truth")
				else:
					outRow.append("Untruth")

				if somaticNess == "NOT_SOMATIC":
					outRow.append("Yes")
					vafMatrixGerm.append(outRow)
				else:
					outRow.append("No")
					vafMatrixSomatic.append(outRow)

			continue

		# back to the normal, usual situation where there is only one alternate allele...
		outRow.append(chrom + ":" + position) 
		
		# get variant allele frequencies - first need to get FORMAT information
		formatInfo = lineList[header.index("FORMAT")].split(":")

		aoIndex = formatInfo.index("AO")
		roIndex = formatInfo.index("RO")

		for sampleIndex in range((formatIndex + 1), len(lineList)):
			sampleRecord = lineList[sampleIndex].split(":")

			AO = sampleRecord[aoIndex]
			RO = sampleRecord[roIndex]

			VAF = ""

			if AO == "" or RO == "" or AO == "." or RO == ".":
				VAF = "NA"
			elif AO == "0" and RO == "0":	
				VAF = "NA"
			else:
				VAF = float(AO) / (float(AO) + float(RO))

			outRow.append(str(VAF))

		# add ref and alt alleles to outRow
		outRow.append(ref)
		outRow.append(alt)
		outRow.append(snpEffInfo)
		outRow.append(mutationType)
		outRow.append(mutationImpact)
		outRow.append(mutationGene)		

		if mutationProteinChange == "":
			outRow.append("None")
		else:
			outRow.append(mutationProteinChange)

		if mutationGene in sangerList:
			outRow.append("Truth")
		else:
			outRow.append("Untruth")

		if somaticNess == "NOT_SOMATIC":
			outRow.append("Yes")
			vafMatrixGerm.append(outRow)
		else:
			outRow.append("No")
			vafMatrixSomatic.append(outRow)

inputFile.close()

# sort matrix by variant impact (3rd-to-last column)
vafMatrixSomatic = sorted(vafMatrixSomatic, key = itemgetter(outHeader.index("Effect"), outHeader.index("SangerGene?")))
vafMatrixGerm = sorted(vafMatrixGerm, key = itemgetter(outHeader.index("Effect"), outHeader.index("SangerGene?")))

# output to file
outputFileSomatic = open(outFileSomatic, "w")

for line in vafMatrixSomatic:
	outputFileSomatic.write("\t".join(line) + "\n")

outputFileSomatic.close()


outputFileGerm = open(outFileGerm, "w")

for line in vafMatrixGerm:
	outputFileGerm.write("\t".join(line) + "\n")

outputFileGerm.close()
