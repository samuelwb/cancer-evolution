## SNV/indel identification and subclone analysis

This folder contains scripts for identifying somatic and germline SNVs and indels
from aligned BAM files (generated using BWA or other approaches; please see "Alignment"
folder in this repository). These scripts also identify and plot mutation clusters
for longitudinal subclone analysis.

Required input includes a germline BAM file for the patient plus 2 or more longitudinal
cancer sample BAMs from the same patient. The inputs for each script are indicated in 
comments at the top of the script itself. Please look at these comments for each script to 
get further details on execution.

The scripts should be executed in the following order:

* **a_FreeBayesChromParallelExecute.sh**: Runs FreeBayes on your files to identify SNVs and
	indels. The output is one VCF file for each chromosome. This VCF file must be combined 
	using an appropriate tool (such as VCFtools) prior to the next step.

* **b_AnnotationSnpEff.sh**: Adds variant annotation to the VCF to indicate variant effects,
	genes affected by the variant, etc.

* **c_QualFilter.sh**: Filters out low-quality variants.

* **d_MarkSomatic.sh**: Marks each variant as somatic or not somatic.

* **e_RemoveNotSomatic.sh**: Removes non-somatic (germline) variants to obtain a somatic-only
	VCF. This step is not necessary but is useful if you would like to manually
	look at only the somatic variants in the VCF.

* **f_MakeVafMatrixFromVcf.sh**: Using VCF output from step d, determines the variant allele
	frequency (VAF) of each mutation for each sample, and outputs this in a tab-delimited
	text file (one file for somatic mutations, one for germline).

* **g_AdjustVaf.sh**: Adjusts somatic VAFs from previous step based on normal contamination (determined
	by copy number analysis or other methods).

* **h_AddCopyNum.sh**: Adds the copy number status (neutral or altered) to a new column in the VAF
	txt file (see "CopyNumberWgs" folder in this repository for how to generate copy number 
	data required for this step). This allows later VAF heatmap plotting with or without 
	copy-altered mutations included.

* **i_ChangePattern.sh**: Groups somatic mutations into mutation clusters based on their pattern of
	change.

* **j_SortVafs.sh**: Sorts VAFs based on change pattern (mutation cluster) for clearer plotting of
	mutation VAFs in heatmaps.

* **k_VafHeatmapAndSubcloneAnalysis.R**: Generates VAF heatmaps based on output of previous step to view
	mutation evolution. This is best run in RStudio.

* **l_MutationSignaturesGermFilter.sh**: Determines trinucleotide mutation context for each SNV and
	summarizes the number of SNVs with each mutation context in each change pattern (mutation cluster).

