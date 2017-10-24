# Author: Stephen Piccolo

library(DNAcopy)

inCopyNumberFilePath = commandArgs()[7]
outCopyNumberFilePath = commandArgs()[8]

set.seed(1)
cn = read.table(inCopyNumberFilePath, header=T, stringsAsFactors=F, sep="\t", row.names=NULL, check.names=F)
CNA.object = CNA(genomdat = cn[,7], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
CNA.smoothed = smooth.CNA(CNA.object)
segs = segment(CNA.smoothed, verbose=0, min.width=2)

write.table(segs$out[,2:6], file=outCopyNumberFilePath, row.names=F, col.names=F, quote=F, sep="\t")

