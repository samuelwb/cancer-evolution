#RNA-Seq expression analysis

This folder contains scripts for determining gene expression from RNA-Seq fastq files, including
both bulk RNA-Seq and scRNA-Seq. These scripts are compatible with either paired- or single-end data. 
The output is files containing TPM, RPKM, and feature counts expression data.

The RSubread.sh script will execute RSubread. It calls the R script AlignRnaSeq.R which is also
found in this folder and performs alignment and transcript counting. Please see the comments
near the top of RSubread.sh for details on what inputs are required to this script.

