# Author: Stephen Piccolo

library(Rsubread)
library(limma)
library(edgeR)
library(tools)
options(digits=2)

referenceGenomeFastaFilePath = commandArgs()[7]
inFilePath1 = commandArgs()[8]
inFilePath2 = commandArgs()[9] # NULL for single-end analyses or when a BAM file has been specified
gtfFilePath = commandArgs()[10]
outFilePrefix = commandArgs()[11]

memory = 16000
nthreads = 24

input_format = "gzFASTQ"
if (file_ext(inFilePath1) == "bam")
  input_format = "BAM"

outBamFilePath = paste(outFilePrefix, "bam", sep=".")
outCountsFilePath = paste(outFilePrefix, "featureCounts", sep=".")
outRpkmFilePath = paste(outFilePrefix, "rpkm", sep=".")
outRpkmLogFilePath = paste(outFilePrefix, "rpkmlog", sep=".")
outTpmFilePath = paste(outFilePrefix, "tpm", sep=".")
outTpmLogFilePath = paste(outFilePrefix, "tpmlog", sep=".")
outStatsFilePath = paste(outFilePrefix, "stats", sep=".")

referenceGenomeIndexFilePrefix = paste(referenceGenomeFastaFilePath, "__reference_index", sep="")

if (!file.exists(paste(referenceGenomeIndexFilePrefix, ".reads", sep="")))
  buildindex(basename=referenceGenomeIndexFilePrefix, reference=referenceGenomeFastaFilePath, memory=memory)

if (inFilePath2 == "NULL")
  inFilePath2 = NULL

if (!file.exists(outBamFilePath))
  align(index=referenceGenomeIndexFilePrefix, readfile1=inFilePath1, readfile2=inFilePath2, output_file=outBamFilePath, nthreads=nthreads, input_format=input_format, tieBreakHamming=TRUE, unique=TRUE, indels=5)

fCountsList = featureCounts(outBamFilePath, annot.ext=gtfFilePath, isGTFAnnotationFile=TRUE, nthreads=nthreads, isPairedEnd=!is.null(inFilePath2))
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
rpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = rpkm / sum(rpkm) * 10^6

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts)
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(cbind(fCountsList$annotation[,1], rpkm), outRpkmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], log2(rpkm + 1)), outRpkmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

write.table(cbind(fCountsList$annotation[,1], tpm), outTpmFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
write.table(cbind(fCountsList$annotation[,1], log2(tpm + 1)), outTpmLogFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

unlink(outBamFilePath)
unlink(paste(outBamFilePath, ".indel", sep=""))
