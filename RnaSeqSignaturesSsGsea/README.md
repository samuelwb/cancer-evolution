## Signature analysis of bulk and single-cell RNA-Seq data

This folder contains scripts for running pathway analysis on single-cell or bulk RNA-Seq
data using ssGSEA. The ssGSEA.sh script calls the ssGSEA.r script. The user can modify
the inputs to the ssGSEA.sh script to run ssGSEA on their particular gene expression data, 
which should be input in the form of a tab-delimited text file with samples in columns
and genes in rows, and gene expression in TPM or FPKM format.

The user can also run the particular set of signatures they are interested in by downloading
a gmt file representing these signatures from http://software.broadinstitute.org/gsea/msigdb.

Please see ssGSEA.sh for details on file inputs to this script.

