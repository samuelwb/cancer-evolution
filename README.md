## About this repository

This repository contains Python, R and other code used in the paper "Combating subclonal evolution of resistant cancer phenotypes" (https://www.ncbi.nlm.nih.gov/pubmed/29093439). In this study, we analyzed the evolution of four ER+ breast cancers over 2-15 years using whole-genome sequencing (WGS), whole-exome sequencing (WES), bulk RNA-Seq, and single-cell RNA-Seq (scRNA-Seq) of longitudinal patient samples. For questions, please contact me (Samuel Brady) at samuelwarrenbrady@gmail.com or Samuel.Brady@utah.edu.

Here you can find our code for analyzing each of these data types, divided into the following subfolders, each of which has its own README:

1.  **Alignment**: Code for aligning DNA fastq reads to the human genome to create an aligned BAM file, providing the basis for variant identification approaches described in other folders.
2.  **CopyNumberWgs**: Code for identifying CNAs from WGS data.
3.  **RnaSeqExpression:** Code for processing bulk RNA-Seq and scRNA-Seq fastq files to obtain gene-level expression value (including RPKM, TPM, and feature counts).
4.  **RnaSeqFusions**: Code for identifying RNA fusion events based on RNA-Seq fastq input.
5.  **RnaSeqSignaturesAssign**: Code for determining gene expression signatures (our custom gene signatures) that change over time in each patient, based on bulk RNA-Seq or scRNA-Seq TPM values obtained with code from the "RnaSeqExpression" folder.
6.  **RnaSeqSignaturesSsGsea**: Code for determining gene expression signatures (public gene signatures from the Molecular Signatures Database) that change over time in each patient, based on scRNA-Seq TPM values obtained with code from the "RnaSeqExpression" folder.
7.  **RnaSeqSingleCellCopyNumber**: Code for inferring CNAs from scRNA-Seq (TPM) data.
8.  **SnvWgs**: Code for identifying SNVs and indels from WGS and for then identifying mutation clusters in longitudinal samples for subclone analysis. Also code for identifying mutation context and mutation signatures of SNVs.
9.  **StructuralWgs**: Code for identifying structural variants (translocations, large deletions, inversions, and duplications) from WGS data.

Please click the subfolder of interest to view the code. Each subfolder contains a README file describing how to run the code in the folder.

