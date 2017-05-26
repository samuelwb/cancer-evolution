# ASSIGN Analysis Scripts

To recreate the ASSIGN analysis performed on the patient data, you can use the
scripts provided here. Please contact [David Jenkins](mailto:dfj@bu.edu) with
any questions or issues.

>Before running the scripts, you must modify the input and output locations
>in both scripts

## R Packages

Before running the scripts, make sure the following R packages are installed:

__CRAN Packages__:

1. [devtools](https://cran.r-project.org/web/packages/devtools/index.html)
2. [data.table](https://cran.r-project.org/web/packages/data.table/index.html)

__Development Packages__:

Development versions of _SVA_ and _ASSIGN_ were used for the analysis in this
paper. Please install the development versions of these packages, available
on GitHub:

```
library(devtools)
install_github("wevanjohnson/sva-devel", ref="reference_batch")
install_github("dfjenkins3/ASSIGN", ref="even_signature")
```

## Running the Scripts

The script ```ASSIGN_merge_and_combat.R``` will load the data, run batch
correction, and save an rsession that can be loaded by
```ASSIGN_run_predictions.R```. Before running either script, modify the file
paths in the _Input Files_ and _Output Files_ section of the scripts to
locations on your system. Once the file paths have been set, you can run the
analysis with the following commands:

```
Rscript ASSIGN_merge_and_combat.R
Rscript ASSIGN_run_predictions.R 1
Rscript ASSIGN_run_predictions.R 2
Rscript ASSIGN_run_predictions.R 3
Rscript ASSIGN_run_predictions.R 4
Rscript ASSIGN_run_predictions.R 5
Rscript ASSIGN_run_predictions.R 6
Rscript ASSIGN_run_predictions.R 7
Rscript ASSIGN_run_predictions.R 8
Rscript ASSIGN_run_predictions.R 9
Rscript ASSIGN_run_predictions.R 10
```

The pathway prediction values for each ASSIGN run can be found in the
```pathway_activity_testset.csv``` file in each pathway's subdirectory.

## Other Files

* ```Key_ASSIGN_functions.R``` - Internal R functions used by
```ASSIGN_merge_and_combat.R``` and ```ASSIGN_run_predictions.R```.
* ```gene_lists.rda``` - Rdata object containing the gene lists used
for the pathways. This file is loaded by ```ASSIGN_run_predictions.R```.
