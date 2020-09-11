library(cmapR)
library(mmalign)
library(dplyr)
library(limma)

# create signature groups


lincs.indir <- ".";

source(file.path(lincs.indir, "common.R"));

# use only landmark genes
pheno <- qread(file.path(lincs.indir, "lincs-sig-info.rds"));


# control groups are not split
controls <- lapply(control.types, function(type) which(pheno$pert_type == type));
names(controls) <- control.types;

pert.types <- levels(pheno$pert_type);
treat.types <- setdiff(pert.types, control.types);

# treatment groups are split by perturbogen gen
treats <- lapply(treat.types,
	function(treat.type) {
		pert.names <- unique(pheno$pert_iname[pheno$pert_type == treat.type]);

		idx <- lapply(pert.names, function(name) which(pheno$pert_iname == name));
		names(idx) <- paste0(treat.type, ",", pert.names);

		idx
	}
);
names(treats) <- treat.types;

groups <- c(list(ctl=controls), treats);

qwrite(groups, "lincs-sig-groups.rds");

