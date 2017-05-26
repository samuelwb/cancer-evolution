# Author: Samuel Brady

# This R script generates a heatmap of VAFs for this patient using the output
# of the previous step (step j). This enables visualization of variant evolution.
# This script requires the ComplexHeatmap and pvclust packages in R and is
# ideally run in RStudio.

# Please modify lines 18 and 20 to indicate your working directory and input file
# (VAF txt file from step j). Please modify line 61 to change the name of your output
# file indicating the VAFs of maximum density for each change pattern (mutation cluster).

library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(MASS)
library(pvclust)

setwd("dir_with_input_file")

vafMatrix <- as.matrix(read.table("output_VAF_matrix_from_step_j", header=TRUE, stringsAsFactors=FALSE))
codingNess <- which(colnames(vafMatrix)=="Ref")
copyIndex <- which(colnames(vafMatrix)=="Log2Copy") # if you want to filter by copy-2 mutations
changePattIndex <- which(colnames(vafMatrix)=="SamplesWithMutation")

vafMatrix[vafMatrix=="Unknown"] <- NA
vafMatrixCut <- vafMatrix[complete.cases(vafMatrix),]
vafMatrixCut <- vafMatrixCut[vafMatrixCut[,copyIndex] == "Neutral", ] # use only if you want to filter for copy-2 mutations
vafMatrixCut <- vafMatrixCut[vafMatrixCut[,2] < 0.001, ] # use only if you want to do more strict germline variant filtering

write.table(vafMatrixCut[,2:(codingNess-1)], "TempSummaryNumeric.txt", quote=FALSE, sep="\t", col.names=NA)
vafMatrixNum <- as.matrix(read.table("TempSummaryNumeric.txt", header=TRUE, stringsAsFactors=FALSE))

Heatmap(vafMatrixNum, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(0,0.5,1), c("cornsilk", "blue", "darkblue")), show_row_names = FALSE) #, split=as.numeric(inMatrix[,ncol(inMatrix)]), 
# Use this section to determine the VAF peaks for each change pattern
# get a list of change patterns
changePatternList <- unique(vafMatrixCut[,changePattIndex])
sampleIndices <- c(2:(codingNess - 1))

chgPttSumm <- data.frame(matrix(NA, nrow = length(changePatternList), ncol = codingNess - 1))

rowCounter <- 1

for (thisChangePattern in changePatternList)
{
  outRow <- c(thisChangePattern)
   
  for (sampleIndex in sampleIndices)
  {
    d <- density(as.numeric(vafMatrixCut[vafMatrixCut[,changePattIndex] == thisChangePattern,][,sampleIndex])); xVal <- d$x[which.max(d$y)]; plot(d); abline(v=xVal); xVal
    outRow <- c(outRow, xVal)
  }
  
  cat("\n\n")
  cat(outRow)
  chgPttSumm[rowCounter,] <- outRow
  rowCounter <- rowCounter + 1
}

colnames(chgPttSumm)[1] <- "ChangePattern"

write.table(chgPttSumm, "MySamples_ChangePatternSummary.txt", quote=FALSE, sep="\t", col.names=NA)

# investigate specific change patterns in greater detail
d <- density(as.numeric(vafMatrixCut[vafMatrixCut[,changePattIndex] == "Samples:1-3-4",][,4])); xVal <- d$x[which.max(d$y)]; plot(d, xlim=c(0,1)); abline(v=xVal); xVal
peakx <- d$x[which(diff(sign(diff(d$y)))==-2)]; peakx

