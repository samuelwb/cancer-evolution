#!/bin/bash

# Authors: Stephen Piccolo and Samuel Brady

# This bash script runs ssGSEA pathway analysis on the input (TPM txt) RNA-Seq expression file

geneExpressionFile=		# txt file with gene expression (gene symbols in rows, samples in columns) in transcripts per million (TPM)
c2="c2.all.v5.2.symbols.gmt"	# download from http://software.broadinstitute.org/gsea/msigdb
outputFile=			# output file with ssGSEA enrichment scores for all gene sets interrogated

Rscript ssGSEA.r $geneExpressionFile $c2 $outputFile
