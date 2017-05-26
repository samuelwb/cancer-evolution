## Inferring copy number from single-cell RNA-Seq data

This folder contains scripts for determining copy number from single-cell RNA-Seq data based
on the algorithm published by Tirosh et al. (Science 352:6282, 2016).

This analysis requires single-cell RNA-Seq data from diploid normal cells as a 2-copy control. We have 
provided human mammary epithelial cell single-cell RNA-Seq in this folder ("HMEC_Combined.tpm"),
which is most suitable for breast cancer data, but it may be better for the user to find
bulk or single-cell RNA-Seq data corresponding to the normal tissue of his/her cancer as 
a better control if this file does not yield suitable results.

The user's input single-cell RNA-Seq data should be prepared as a tab-delimited text file with "\n"
characters (not "\r") ending each line, with genes in rows and samples (single cells) in columns,
and data in TPM (transcripts per million) format.

The scripts should be executed in the following order (please see comments at the top of each bash
script for details about how to execute the script):

* **a_CombineWithHmec.sh**: Combines the user's single-cell RNA-Seq gene expression txt file with
	normal cell (such as HMEC) expression data, which is added on as additional columns,
	for normalization.

* **b_SingleCellCna.sh**: Infers the copy number of each single-cell using 101-gene window expression
	averages.

* **c_RemoveNormal.sh**: Removes normal cells used for normalization from the previous step's output.

* **d_PostNormalization.sh**: After the previous step ("c_RemoveNormal.sh") the user should plot the
	inferred copy number of each cell in a heatmap showing each chromosome, and identify
	cells lacking copy number alterations (normal cells) and cells with copy number alterations
	matching bulk DNA-Seq inferred copy number (cancer cells). Then make a text file patterned
	after "cell_identity.txt" found in this folder to indicate whether each cell is normal or cancer.
	This file is used as input to this step ("d_PostNormalization.sh") to perform additional
	normalization to normal cells *within the user's experiment* which will reduce noise associated
	with using normal-cell gene expression data from a different experiment (i.e. HMECs). After
	this step the user will have inferred copy for each single-cell and can plot the copy number
	of each single cell's chromosomes using a heatmap or other desired approaches.


