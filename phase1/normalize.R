library(io)
library(cmapR)
library(dplyr)
library(ggplot2)
library(reshape2)

lincs.indir <- "~/data/cmap/phase1";

source(file.path(lincs.indir, "common.R"))

gct <- file.path(lincs.indir, "lincs-sig_lm_n473647x978.gctx");

sig <- qread("lincs-sig-info.rds");
gene <- qread("lincs-gene-info_lm.rds");
cell <- qread("lincs-cell-info.rds");

# sanity check

rid <- read.gctx.ids(gct);
cid <- read.gctx.ids(gct, dimension="column");

n <- length(cid);

stopifnot(as.integer(rid) == gene$pr_gene_id);
stopifnot(cid == sig$sig_id);


# read all consensus control types

sig.rc <- filter(sig, pert_type %in% c.control.types);

#rc <- parse.gctx(gct, cid = sig.rc$sig_id)@mat;

stopifnot(sig.rc$sig_id == colnames(rc))

# @param x   data matrix of samples to be normalized
# @param sig sample information of data matrix
# @param rc  data matrix of consensus reference samples
normalize_gctx <- function(x, sig, rc) {
	ids.x <- colnames(x);

	sig.x <- sig[match(ids.x, sig$sig_id), ];
			
	# FIXME slow!
	matched <- lapply(1:nrow(sig.x),
		function(i) {
			best_matched_control(sig.x[i, ], sig)
		}
	);
	
	mapply(
		function(treatment, controls) {
			x[, treatment] - rowMeans(rc[, controls[controls %in% colnames(rc)], drop=FALSE])
		},
		ids.x,
		matched
	)
}


# test on one chunk

#x <- parse.gctx(gct, cid = 1:100, matrix_only=TRUE)@mat;
#xn <- normalize_gctx(x, sig, rc);
#plot(x, xn, pch=".")
#cor(c(x), c(xn))

xns <- apply_gctx(gct, normalize_gctx, n, sig=sig, rc=rc);

xn <- do.call(cbind, xns);

ds <- methods::new("GCT", mat = xn, rid = NULL, cid = NULL, 
	set_annot_rownames = FALSE, matrix_only=TRUE);

out.fname <- "lincs-sig_lm_norm";
write.gctx(ds, out.fname, matrix_only=TRUE);

