library(io)
library(cmapR)
library(dplyr)
library(ggplot2)
library(reshape2)

lincs.indir <- "~/data/cmap/phase1";

gct <- file.path(lincs.indir, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx");

sig <- qread("lincs-sig-info.rds");
gene <- qread("lincs-gene-info.rds");
cell <- qread("lincs-cell-info.rds");

# sanity check

rid <- read.gctx.ids(gct);
cid <- read.gctx.ids(gct, dimension="column");

stopifnot(as.integer(rid) == gene$pr_gene_id);
stopifnot(cid == sig$sig_id);

#

dmso <- filter(sig, pert_iname == "DMSO");
dmso.cns <- filter(sig, pert_iname == "DMSO", pert_type == "ctl_vehicle.cns");
dox <- filter(sig, pert_iname == "doxorubicin");

x <- parse.gctx(gct, cid = dox$sig_id, matrix_only=TRUE)@mat;
r <- parse.gctx(gct, cid = dmso$sig_id, matrix_only=TRUE)@mat;
rc <- parse.gctx(gct, cid = dmso.cns$sig_id, matrix_only=TRUE)@mat;


qdraw(
	{
		plot(x[,4], type="l")
	},
	file = "lincs-dox4.pdf"
);

plot(x[,1], type="l")
plot(x[,2], type="l")
plot(x[,3], type="l")
plot(x[,5], type="l")

plot(r[,1], type="l")
lines(gene$pr_is_bing, col="red")
plot(r[,2], type="l")
plot(r[,3], type="l")
plot(r[,4], type="l")
plot(r[,5], type="l")

blocks <- c(rep(0, 1000), rep(1, 6000), rep(2, nrow(x) - 7000));

qdraw(
	{
		par(mfrow=c(3,1));
		hist(x[blocks == 0, ], breaks=1000);
		hist(x[blocks == 1, ], breaks=1000);
		hist(x[blocks == 2, ], breaks=1000);
	},
	width = 4, height = 8,
	file = "lincs-blocks_dox_hist.pdf"
)

qdraw(
	{
		par(mfrow=c(3,1));
		hist(r[blocks == 0, ], breaks=1000);
		hist(r[blocks == 1, ], breaks=1000);
		hist(r[blocks == 2, ], breaks=1000);
	},
	width = 4, height = 8,
	file = "lincs-blocks_dmso_hist.pdf"
)

qdraw(
	{
		par(mfrow=c(3,1));
		hist(rc[blocks == 0, ], breaks=500, fill="black");
		hist(rc[blocks == 1, ], breaks=500, fill="black");
		hist(rc[blocks == 2, ], breaks=500, fill="black");
	},
	width = 4, height = 8,
	file = "lincs-blocks_dmso-cns_hist.pdf"
)

d <- melt(x, varnames=c("gene_id", "sig_id")) %>%
	left_join(select(sig, -distil_id));

ggplot(d, aes(x=value, fill=cell_id)) +
	geom_density(alpha=0.2) + theme_bw() +
	facet_grid(pert_idose ~ .)

