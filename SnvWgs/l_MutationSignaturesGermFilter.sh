#!/bin/bash

# Author: Samuel Brady

# This bash script can be executed on the output file from step j to add a column to
# your VAF txt file that indicates the trinucleotide mutation context of each SNV. It
# also summarizes the number of mutations in each mutation cluster (each change pattern)
# with each mutation context (96 possible mutation contexts).

# This data can be used as input to deconstructSigs to perform mutation signature analysis.
# Please see documentation for deconstructSigs or contact Samuel Brady (Samuel.Brady@utah.edu)
# for more information on how to interface this output with deconstructSigs.

inVafMatrix=											# input VAF matrix, which was output from step j
chromFastaFolder=/data/Genomes/iGenomes/Homo_sapiens/Ensembl/GRCh37/Sequence/Chromosomes	# folder containing fasta sequences for each chromosome in the form "1.fa", "2.fa" for chromosomes 1 and 2, etc.
germVafFilter=0.001										# exclude variants germline VAF above this value
outVafMatrix=											# output file with VAF matrix with mutation context added as a new colum
mutSigSummaryFile=79NTL6_MutSigSummaryG.txt							# output file summarizing number of SNVs with each mutation context in each change pattern (mutation cluster)

python MutationSignaturesGermFilter.py $inVafMatrix $chromFastaFolder $germVafFilter $outVafMatrix $mutSigSummaryFile

