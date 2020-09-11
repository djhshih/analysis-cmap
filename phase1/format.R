library(io)
library(cmapR)

lincs.indir <- "~/data/cmap/phase1";

gct <- file.path(lincs.indir, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx");

pheno <- qread(file.path(lincs.indir, "GSE92742_Broad_LINCS_sig_info.txt"), type="tsv", stringsAsFactors=FALSE, na.strings="-666");

fannot <- qread(file.path(lincs.indir, "GSE92742_Broad_LINCS_gene_info.txt"), type="tsv", stringsAsFactors=FALSE, na.strings="-666");

cell <- qread(file.path(lincs.indir, "GSE92742_Broad_LINCS_cell_info.txt"), type="tsv",
	stringsAsFactors=TRUE, na.strings="-666");

# Reformat data ##############################################################

# rearrange fannot to match data

rid <- as.integer(read.gctx.ids(gct));
cid <- read.gctx.ids(gct, dimension="column");

stopifnot( !any(is.na(rid)) )

length(rid)
length(cid)

idx <- match(rid, fannot$pr_gene_id);
table(is.na(idx))
fannot <- fannot[idx, ];

stopifnot( all(fannot$pr_gene_id == rid) )

# rearrange pheno to match data

idx <- match(cid, pheno$sig_id);
table(is.na(idx))
pheno <- pheno[idx, ];

stopifnot( all(pheno$sig_id == cid) )


# remove redundant columns

stopifnot( with(pheno, all(paste(pert_time, pert_time_unit) == pert_itime)))
pheno$pert_itime <- NULL;

# annotate missing data properly
idx <- pheno$pert_dose == "-666.0|-666.000000" | pheno$pert_dose == "-666.000000";
pheno$pert_dose[idx] <- NA;
pheno$pert_idose[idx] <- NA;

# fix incorrect entries
idx <- pheno$pert_dose == "300.0|300.000000"
pheno$pert_dose[idx] <- 300.0;
pheno$pert_idose[idx] <- "300 ng";

pheno$pert_dose <- as.numeric(pheno$pert_dose);
pheno$pert_dose[pheno$pert_dose < 0] <- NA;

# paste(pert_dose, pert_dose_unit) != pert_idose
# presumably because pert_dose contains the actual dose (due to purity?)
# and pert_idose indicates the target dose?

# convert selected fields to factors
pheno$cell_id <- factor(pheno$cell_id);
pheno$pert_dose_unit <- factor(pheno$pert_dose_unit);
pheno$pert_time_unit <- factor(pheno$pert_time_unit);
pheno$pert_type <- factor(pheno$pert_type);

qwrite(pheno, "lincs-sig-info.tsv");
qwrite(pheno, "lincs-sig-info.rds");

qwrite(fannot, "lincs-gene-info.tsv");
qwrite(fannot, "lincs-gene-info.rds");

qwrite(cell, "lincs-cell-info.tsv");
qwrite(cell, "lincs-cell-info.rds");


# subset only landmark genes
idx <- which(fannot$pr_is_lm == 1);

fannot.lm <- fannot[idx, ];
ds.lm <- parse.gctx(gct, rid=idx, matrix_only=TRUE);

qwrite(fannot.lm, "lincs-gene-info_lm.tsv");
qwrite(fannot.lm, "lincs-gene-info_lm.rds");

write.gctx(ds.lm, "lincs-sig_lm", matrix_only=TRUE);

