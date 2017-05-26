#!/bin/bash

# Author: Samuel Brady

# This script executes TopHat-Post, the second step in identifying fusion
# transcripts from RNA-Seq data. TopHat is very sensitive to correct
# folder structure. This script should be executed in the folder containing
# the "tophat_[MySampleName]" output folder from the previous step. This
# script will automatically identify any folders named in this way in the
# current folder and will perform analysis. The output data will then
# be placed into a folder called "tophatfusion_out."

# User-defined inputs
fusionDatabasePath=     # path to database files included in TopHat installation (i.e. "/data/Software/tophat-2.0.6.Linux_x86_64/TopHat_Fusion_Database_Files/hg19")

# Transferring blast files included in the TopHat installation to the current directory may be necessary for successful implementation
echo "Transferring blast files to this directory; I will move them back afterwards."
mv /path/to/tophat-2.0.6.Linux_x86_64/top_dir/blast ./

echo "---Beginning TopHat Post---"
/path/to/tophat-2.0.6.Linux_x86_64/tophat-fusion-post -p 8 --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 $fusionDatabasePath

echo "Transferring blast files back to original location"
mv blast /path/to/tophat-2.0.6.Linux_x86_64/
