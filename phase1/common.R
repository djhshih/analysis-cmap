library(cmapR)

# find best matching control
best_matched_control <- function(case.sig.entry, control.sig.db) {
	type <- case.sig.entry$pert_type;

	control.cell <- as.character(case.sig.entry$cell_id);
	control.time <- as.character(case.sig.entry$pert_time);
	control.dose <- as.character(case.sig.entry$pert_idose);

	if (type %in% c("trt_cp", "ctl_vehicle", "ctl_vehicle.cns")) {
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
	} else if (type %in% c("trt_sh", "trt_sh.cgs", "trt_sh.css")) {
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

# apply a function to column chunks from a gctx file
# @param n  number of columns in gctx file
apply_gctx <- function(fname, f, n, chunk.size=10000, ...) {
	require(cmapR)

	bstarts <- seq(1, n, by=chunk.size);

	lapply(bstarts,
		function(bstart) {
			message("block start ", bstart)
			istart <- bstart;
			iend <- min(bstart + chunk.size - 1, n);

			mat <- parse.gctx(fname, cid=istart:iend, matrix_only=TRUE)@mat;
			
			f(mat, ...)
		}
	)
}

c.control.types <- c("ctl_vehicle.cns", "ctl_untrt.cns", "ctl_vector.cns", "trt_sh.css");
control.types <- c("ctl_vehicle", "ctl_vehicle.cns", "ctl_untrt", "ctl_untrt.cns", "ctl_vector", "ctl_vector.cns", "trt_sh.css");

