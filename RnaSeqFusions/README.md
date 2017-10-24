## RNA-Seq fusion transcript identification

This folder contains code for identifying fusion transcripts
from RNA-Seq data. It requires installation of TopHat. (We used
TopHat v2.0.6.) The fusion search is performed in two steps (two
bash scripts), the latter of which is called TopHat-Post. Each of
these bash scripts has a section at the top of the script indicating
the required inputs.

The scripts should be executed in the following order:

* **a_TopHat2.sh**: Performs preliminary steps required for fusion identification.

* **b_TopHat2_Post.sh**: Analyzes output from previous step and generates html
	output showing fusion transcripts identified. See comments at the
	top of the script for further information.

