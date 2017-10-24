# Authors: Stephen Piccolo and Samuel Brady

# This script performs ssGSEA on RNA seq data (matrix)
source("http://bioconductor.org/biocLite.R")
library(GSEABase)
library(GSVAdata)
data(c2BroadSets)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)

cat("\n\n***Successfully loaded required libraries***\n\n")

# Get command-line arguments
args <- commandArgs(TRUE)
inputfilename <- args[1]
gmtFile <- args[2]
outputFileName <- args[3]

# Read in RNA-Seq matrix
RNAseqMatrix <- as.matrix(read.table(inputfilename, header=TRUE, stringsAsFactors=FALSE, row.names=1))

# Within GSEABase, function converts gmt files to an object that has those gene sets (convert gene set collection class, because
# the c2BroadSets gene set collection class uses gene IDs instead of gene symbols); use the getGMT
c2set <- getGmt(gmtFile)

print(c2set)

# Run GSVA in ssGSEA format and output to file
gsvaData <- gsva(RNAseqMatrix, c2set, min.sz=10, max.sz=500, verbose=TRUE, rnaseq=TRUE, method="ssgsea")
write.table(gsvaData, outputFileName, col.names=FALSE, quote=FALSE, sep="\t")

