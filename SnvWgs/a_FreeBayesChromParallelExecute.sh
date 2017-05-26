#!/bin/bash

# Author: Samuel Brady

# This bash script executes FreeBayes on aligned bam files to identify
# SNVs and indels (see "Alignment" folder in this repository). 
# Please download FreeBayes from https://github.com/ekg/freebayes.
# This bash script uses the "parallel" commmand to run SNV and indel 
# calling on your aligned bam files on each chromosome in parallel. 

# It executes the commands found in the txt file
# "FreeBayesChromParallel.txt" file in parallel; each line of that file
# performs variant calling on one chromosome for your sample list. It has
# 25 lines, so running this command with the 25 argument at the end will
# execute all chromosomes at the same time, but reducing the number will
# work also (chromosomes will be run in batches in that case).

# You will need to change the arguments in FreeBayesChromParallel.txt,
# which is included in this folder, in the following ways:
#
# 1. Make the "-f" argument point to your genome index. 
# 2. Make the "-L" argument point to a txt file containing a list of
#	bams on which FreeBayes should be run. An example file ("FileList.txt")
#	is included in this folder. This should include the germline sample first,
#	followed by cancer samples for this patient in chronological order (for
#	longitudinal analysis of cancer evolution, to which our analysis is geared). 
# 3. Modify the "-v" argument, if desired, to change the output file names.
# 4. Modify the "-r" arguments, if you are not using reference genome GRCh37, to
#	include the correct chromosome ranges.
#
# It is recommended that each patient's samples are run in separate batches. Thus
# this script would be executed once for each patient.

parallel -a FreeBayesChromParallel.txt --ungroup --max-procs 25

