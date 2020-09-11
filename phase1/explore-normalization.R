library(io)
library(cmapR)
library(dplyr)
library(ggplot2)
library(reshape2)

lincs.indir <- "~/data/cmap/phase1";

gct <- file.path(lincs.indir, "lincs-sig_lm_n473647x978.gctx");

sig <- qread("lincs-sig-info.rds");
gene <- qread("lincs-gene-info_lm.rds");
cell <- qread("lincs-cell-info.rds");

# sanity check

rid <- read.gctx.ids(gct);
cid <- read.gctx.ids(gct, dimension="column");

stopifnot(as.integer(rid) == gene$pr_gene_id);
stopifnot(cid == sig$sig_id);

#

dmso <- filter(sig, pert_iname == "DMSO", pert_type == "ctl_vehicle");
dmso.cns <- filter(sig, pert_iname == "DMSO", pert_type == "ctl_vehicle.cns");
dox <- filter(sig, pert_iname == "doxorubicin");

x <- parse.gctx(gct, cid = dox$sig_id, matrix_only=TRUE)@mat;
r <- parse.gctx(gct, cid = dmso$sig_id, matrix_only=TRUE)@mat;
rc <- parse.gctx(gct, cid = dmso.cns$sig_id, matrix_only=TRUE)@mat;

# this will assume > 16 G of memory and cause system to hang
#g <- parse.gctx(gct, matrix_only=TRUE)@mat;

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

hist(x, breaks=100)
hist(r, breaks=100)
hist(r[2,], breaks=100)
hist(rc, breaks=100)

r1 <- r[1,];
qqnorm(r1)
hist(r1, freq=FALSE, breaks=100)
curve(dnorm(x, mean(r1), sd(r1)), col="red", add=TRUE)
curve(dt((x - mean(r1))/sd(r1), df=10), col="blue", add=TRUE)

r.sd <- apply(r, 1, sd);
hist(r.sd, breaks=100)
summary(r.sd)

plot(r.sd)

r.mean <- apply(r, 1, mean);
hist(r.mean, breaks=100)
summary(r.mean)


d1 <- melt(x, varnames=c("gene_id", "sig_id")) %>%
	left_join(select(sig, -distil_id));

ggplot(d1, aes(x=value, fill=cell_id)) +
	geom_density(alpha=0.2) + theme_bw() +
	facet_grid(pert_idose ~ .)

d2 <- melt(rc, varnames=c("gene_id", "sig_id")) %>%
	left_join(select(sig, -distil_id));

d <- rbind(d1, d2);

ggplot(d, aes(x=value, fill=pert_type)) +
	geom_density(alpha=0.2) + theme_bw()
#	facet_grid(cell_id ~ .)
														

# find best matching control
best_matched_control <- function(case.sig.entry, control.sig.db) {
	type <- case.sig.entry$pert_type;

	control.cell <- as.character(case.sig.entry$cell_id);
	control.time <- as.character(case.sig.entry$pert_time);
	control.dose <- as.character(case.sig.entry$pert_idose);

	if (type %in% c("trt_cp", "trt_vehicle", "trt_vehicle.cns")) {
		control.type <- "ctl_vehicle.cns";
		# dose does not need to match
		# (dose of drug in control is 0; DMSO is delivered at fixed dose)
		control.dose <- NA;
	} else if (type %in% c("trt_lig", "ctl_untrt", "ctl_untrt.cns")) {
		control.type <- "ctl_untrt.cns";
		control.dose <- NA;
	} else if (type %in% c("trt_oe", "trt_oe.mut", "ctl_vector", "ctl_vector.cns")) {
		control.type <- "ctl_vector.cns";
		control.dose <- NA;
	} else if (type == c("trt_sh", "trt_sh.cgs", "trt_sh.css")) {
		# consensus controls sharing same seed sequence
		# NB there will be multiple ones
		control.type <- "trt_sh.css";	
		# dose should match, if possible
		# is this siRNA instead of shRNA?
	} else {
		stop("Unknown signature type: ", type);
	}

	# NA with any value returns TRUE
	are_matched <- function(x, y) {
		if (is.na(x) || is.na(y)) TRUE else x == y;
	}


	m <- filter(control.sig.db,
		are_matched(pert_type, control.type),
		are_matched(cell_id, control.cell),
		are_matched(pert_time, control.time),
		are_matched(pert_idose, control.time),
	)$sig_id;

	# relax the criteria, if necessary

	if (length(m) == 0) {
		m <- filter(control.sig.db,
			are_matched(pert_type, control.type),
			are_matched(cell_id, control.cell),
			are_matched(pert_time, control.time)
		)$sig_id;

		if (length(m) == 0) {
			m <- filter(control.sig.db,
				are_matched(pert_type, control.type),
				are_matched(cell_id, control.cell)
			)$sig_id;

			if (length(m) == 0) {
				m <- filter(control.sig.db,
					are_matched(pert_type, control.type)
				)$sig_id;
			}
		}
	}
	
	m
}


sig.ids <- colnames(x);
sig.dox <- sig[match(colnames(x), sig$sig_id), ];

matched <- lapply(1:nrow(sig.dox),
	function(i) {
		best_matched_control(sig.dox[i, ], sig)
	}
);

xn <- mapply(
	function(treatment, controls) {
		message("treatment ", treatment, "; control ", controls);
		x[, treatment] - rowMeans(rc[, controls[controls %in% colnames(rc)], drop=FALSE])
	},
	sig.ids,
	matched
);

plot(x, xn, pch=".")
cor(c(x), c(xn))

rc.cor <- cor(rc);
diag(rc.cor) <- NA;
hist(rc.cor, breaks=100)
x.rc.cor <- cor(x, rc);
hist(x.rc.cor, breaks=100)

