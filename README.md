# Transcriptomic analysis with the Broad Connective Map

This is a collection of scripts for retrieving (phase 1) CMap data and
formatting it for use with R.

## Installation

Clone this repository. From within the `phase1` directory, do

```
./get.sh
Rscript format.R
Rscript sig-groups.R
```

Only the landmark genes are used because the imputed expression for the other
genes exhibit discernible differences in distributions.

Now, we will have the following data files:

```
lincs-sig_lm_n473647x978.gctx      expression values for landmark genes
lincs-sig-info.rds                 signature annotation
lincs-cell-info.rds                cell line annotation
lincs-gene-info.rds                gene annotation
lincs-sig-groups.rds               categorization of signatures
```

Additionally, the `common.R` script contains functions that are useful for
working with this data.

