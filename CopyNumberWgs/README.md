## Copy number analysis

This folder contains scripts for identifying copy number alterations (CNAs) from
tumor WGS data. It requires an aligned BAM file for each sample as input. (See "Alignment"
folder for code to generate aligned BAM files using BWA in the Speed-Seq suite.)
The inputs for each script are indicated in comments at the top of the script itself.

Two of the scripts (steps f and g) adjust tumor copy number data for normal contamination so 
that absolute copy peaks will be centered at 1, 2, 3, 4, etc. instead of shifted towards
2, which occurs in samples with even minimal normal contamination. This process is described
schematically below. This process aids in determining tumor purity and is essential for
identifying subclonal CNAs (with copy number between integer values).

![picture](img/CopyShift.png | width=100)

The scripts should be executed in the following order:

* **a_CnaPileup.sh**: Generates a pileup based on aligned BAM input. Requires samtools (such as
	version 0.1.19).

* **b_RemoveZeroLines.sh**: Based on pileup from previous step, creates a new pileup with lines
	with 0 copy removed, since such lines confound downstream steps. It calls a
	Python (2.7) script also found in this folder.

* **c_VarScanCopy.sh**: Using previous script output, determines copy number using VarScan. Requires
	VarScan (we used VarScan v2.3.9) and Java.

* **d_CopyCaller.sh**: Using previous script output, calls VarScan to further refine copy number
	determination.

* **e_Segment.sh**: Calls an R script (also found in this folder) to perform segmentation on the
	output of the previous step.

* **f_CenterSegmentedAtZero.sh**: Calls a Python (2.7) script, also found in this folder, to center
	the segmented values around zero. 

* **g_OneThreeFourCent.sh**: Calls a Python (2.7) script, also found in this folder, to convert log2
	fold change from 0 copy values to absolute copy and multiply these scores around the 2-axis
	(diploid) axis to make copy-loss and copy-gain regions near 1, 3, and 4. This is necessary
 	since, due to normal contamination present in nearly all tumor samples, the CNA peaks for
	these copy-loss and copy-gain regions tend to collapse towards 2 and not be right on 1, 
	3, or 4 as they should be. Thus this step essentially adjusts CNA data for normal contamination.

* **h_WindowCna.sh**: Calls a Python (2.7) script, also found in this folder, to take 30-segment
	window averages (or number of segments of user's choosing) of CNA data for noise reduction.

* **i_GetGeneLevelCopyNum.sh**: Optional script that acts on the output of e_Segment.sh to produce
	gene-level segmented CNA data. Calls a Python (2.7) script found in this folder.

