# Tests for the MCMC / Bayesian network step (R/inferCNV_BayesNet.R and the
# BUGS models in inst/).
#
# rjags seeding: base R set.seed() does not reach JAGS, and whether
# run_gibb_sampling() seeds its chains is not assumed here. Every test that
# runs JAGS therefore asserts only properties of the posterior (shape, ranges,
# sums, conjugate expectations with a tolerance sized from the 3000 retained
# draws, state assignment of widely separated data), which hold with or
# without seeded chains.

bn_skip_if_no_jags <- function() {
	testthat::skip_if_not_installed("rjags")
	testthat::skip_if_not_installed("coda")
}

# Fixture: 12 genes x 8 cells. Cells 1-4 are the reference (1 +/- 0.05
# alternating), cells 5-8 are tumor. Region 1 (genes 1-6) sits in cells 5-6 and
# is amplified by the i3 mean delta (the HMM state recorded for it is 1, i.e.
# deliberately wrong); region 2 (genes 7-12) sits in cells 7-8 at the normal
# level (the HMM state recorded is 3, also wrong).
bn_ref_sd <- function() {
	s <- 0.05
	s * sqrt(48 / 47)  # sd of 48 values equal to 1 +/- 0.05 in equal numbers
}

bn_delta <- function() {
	abs(qnorm(0.05, mean = 0, sd = bn_ref_sd()))
}

bn_make_expr <- function() {
	genes <- paste0("g", 1:12)
	cells <- paste0("c", 1:8)
	m <- matrix(1, nrow = 12, ncol = 8, dimnames = list(genes, cells))
	for (i in 1:12) {
		for (j in 1:4) {
			m[i, j] <- 1 + 0.05 * (-1)^(i + j)
		}
	}
	m[1:6, 5:6] <- 1 + bn_delta()
	m
}

bn_make_infercnv <- function(m = bn_make_expr()) {
	gene_order <- data.frame(chr = rep("chr1", nrow(m)),
	                         start = seq_len(nrow(m)),
	                         stop = seq_len(nrow(m)),
	                         row.names = rownames(m))
	methods::new("infercnv",
	             expr.data = m,
	             count.data = m,
	             gene_order = gene_order,
	             reference_grouped_cell_indices = list(normal = 1:4),
	             observation_grouped_cell_indices = list(tumor = 5:8))
}

bn_write_fixture_dir <- function(dir) {
	dir.create(dir, showWarnings = FALSE, recursive = TRUE)
	groups <- data.frame(cell_group_name = c("grpA", "grpA", "grpB", "grpB"),
	                     cell = c("c5", "c6", "c7", "c8"))
	utils::write.table(groups, file.path(dir, "17_HMM_pred.cell_groupings"),
	                   sep = "\t", quote = FALSE, row.names = FALSE)
	pred <- data.frame(cell_group_name = c(rep("grpA", 6), rep("grpB", 6)),
	                   gene_region_name = c(rep("chr1-region_1", 6), rep("chr1-region_2", 6)),
	                   state = c(rep(1, 6), rep(3, 6)),
	                   gene = paste0("g", 1:12))
	utils::write.table(pred, file.path(dir, "17_HMM_pred.pred_cnv_genes.dat"),
	                   sep = "\t", quote = FALSE, row.names = FALSE)
	invisible(dir)
}

bn_args <- function(dir, out_dir, HMM_type = "i3", ...) {
	model <- if (HMM_type == "i3") "BUGS_Mixture_Model_i3" else "BUGS_Mixture_Model"
	a <- list(file_dir = dir,
	          model_file = system.file(model, package = "infercnv"),
	          CORES = 1,
	          out_dir = out_dir,
	          resume_file_token = "",
	          plotingProbs = FALSE,
	          postMcmcMethod = "removeCNV",
	          quietly = TRUE,
	          BayesMaxPNormal = 0,
	          HMM_type = HMM_type,
	          k_obs_groups = 1,
	          cluster_by_groups = FALSE,
	          reassignCNVs = TRUE,
	          diagnostics = FALSE)
	extra <- list(...)
	for (nm in names(extra)) a[[nm]] <- extra[[nm]]
	a
}

# object built through initializeObject() from the fixture files, with the
# i3 mean/precision of the reference set by MeanSD()
bn_make_obj <- function(m = bn_make_expr(), set_mean_sd = TRUE) {
	dir <- tempfile("bn_fixture")
	bn_write_fixture_dir(dir)
	inf <- bn_make_infercnv(m)
	obj <- methods::new("MCMC_inferCNV")
	suppressMessages(capture.output(
		obj <- initializeObject(obj, bn_args(dir, tempfile("bn_out")), inf)))
	if (set_mean_sd) {
		suppressMessages(capture.output(
			obj <- MeanSD(obj, HMM_states = NULL, infercnv_obj = inf)))
	}
	obj
}

# object with hand-set slots; no JAGS and no files
bn_manual_obj <- function(cnv_probs, state, cell_probs = NULL, HMM_type = "i3", ...) {
	n <- length(cnv_probs)
	obj <- methods::new("MCMC_inferCNV")
	a <- list(HMM_type = HMM_type, BayesMaxPNormal = 0.5, quietly = TRUE,
	          out_dir = tempfile("bn_manual_out"), plotingProbs = FALSE,
	          postMcmcMethod = "removeCNV", reassignCNVs = TRUE)
	extra <- list(...)
	for (nm in names(extra)) a[[nm]] <- extra[[nm]]
	obj@args <- a
	obj@cell_gene <- lapply(seq_len(n), function(i) {
		list(cnv_regions = factor(paste0("r", i)), Genes = (2 * i - 1):(2 * i),
		     Cells = (2 * i - 1):(2 * i), State = state[i])
	})
	obj@cnv_probabilities <- cnv_probs
	obj@cell_probabilities <- if (is.null(cell_probs)) as.list(seq_len(n)) else cell_probs
	obj
}

# N draws whose column means are exactly `means` (alternating +/- 0.01)
bn_theta_draws <- function(means, n = 4) {
	p <- matrix(rep(means, each = n), nrow = n)
	p[, 1] <- p[, 1] + c(0.01, -0.01)
	p[, 2] <- p[, 2] - c(0.01, -0.01)
	colnames(p) <- paste0("theta[", seq_along(means), "]")
	p
}

bn_sim_obj <- function(mu, sig, HMM_type = "i3") {
	obj <- methods::new("MCMC_inferCNV")
	obj@mu <- mu
	obj@sig <- sig
	obj@args <- list(quietly = TRUE, HMM_type = HMM_type)
	obj@bugs_model <- system.file(if (HMM_type == "i3") "BUGS_Mixture_Model_i3" else "BUGS_Mixture_Model",
	                              package = "infercnv")
	obj
}

# ---------------------------------------------------------------- accessors

test_that("accessors return the corresponding slots", {
	obj <- methods::new("MCMC_inferCNV")
	obj@cell_gene <- list(list(cnv_regions = "a", Genes = 1:2, Cells = 3L, State = 1))
	obj@bugs_model <- "some/model/file"
	obj@args <- list(HMM_type = "i3", CORES = 2)
	obj@expr.data <- matrix(1:6, 2, 3)
	expect_identical(cellGene(obj), obj@cell_gene)
	expect_identical(modelFile(obj), "some/model/file")
	expect_identical(getArgs(obj), list(HMM_type = "i3", CORES = 2))
	expect_identical(expirData(obj), matrix(1:6, 2, 3))
})

test_that("setBayesMaxPNormal sets the threshold and keeps the other arguments", {
	obj <- methods::new("MCMC_inferCNV")
	obj@args <- list(HMM_type = "i3", BayesMaxPNormal = 0)
	obj2 <- setBayesMaxPNormal(obj, 0.37)
	expect_identical(getArgs(obj2)$BayesMaxPNormal, 0.37)
	expect_identical(getArgs(obj2)$HMM_type, "i3")
	# value semantics: the input object is unchanged
	expect_identical(getArgs(obj)$BayesMaxPNormal, 0)
	# also works when the argument is absent
	obj3 <- methods::new("MCMC_inferCNV")
	obj3@args <- list(HMM_type = "i6")
	expect_identical(getArgs(setBayesMaxPNormal(obj3, 0.9))$BayesMaxPNormal, 0.9)
})

# -------------------------------------------------------------- BUGS models

test_that("the BUGS model files exist and declare the expected number of states", {
	f6 <- system.file("BUGS_Mixture_Model", package = "infercnv")
	f3 <- system.file("BUGS_Mixture_Model_i3", package = "infercnv")
	expect_true(nzchar(f6) && file.exists(f6))
	expect_true(nzchar(f3) && file.exists(f3))
	t6 <- paste(readLines(f6, warn = FALSE), collapse = "\n")
	t3 <- paste(readLines(f3, warn = FALSE), collapse = "\n")
	expect_true(grepl("theta[1:6] ~ ddirich(alpha[])", t6, fixed = TRUE))
	expect_true(grepl("theta[1:3] ~ ddirich(alpha[])", t3, fixed = TRUE))
	expect_true(grepl("mu[6]", t6, fixed = TRUE))
	expect_false(grepl("mu[4]", t3, fixed = TRUE))
	expect_true(grepl("mu[3]", t3, fixed = TRUE))
	expect_true(grepl("for(i in 1:6)", t6, fixed = TRUE))
	expect_true(grepl("for(i in 1:3)", t3, fixed = TRUE))
})

# ---------------------------------------------------------- cnv_prob / cell_prob

test_that("cnv_prob keeps exactly the theta columns in order", {
	m <- cbind(`epsilon[1]` = c(1, 2), `epsilon[2]` = c(3, 1),
	           `theta[1]` = c(0.2, 0.3), `theta[2]` = c(0.3, 0.3), `theta[3]` = c(0.5, 0.4))
	out <- cnv_prob(m)
	expect_equal(dim(out), c(2L, 3L))
	expect_identical(colnames(out), c("theta[1]", "theta[2]", "theta[3]"))
	expect_equal(unname(out), unname(m[, 3:5]))
})

bn_ref_cell_prob <- function(m, K) {
	eps <- m[, grepl("epsilon", colnames(m)), drop = FALSE]
	out <- matrix(0, K, ncol(eps))
	for (j in seq_len(ncol(eps))) {
		for (k in seq_len(K)) out[k, j] <- sum(eps[, j] == k) / nrow(eps)
	}
	out
}

test_that("cell_prob gives per-cell state frequencies (i3), worked by hand", {
	obj <- bn_manual_obj(list(), numeric(0), HMM_type = "i3")
	# 4 draws, 2 cells: cell 1 states 1,1,2,3 ; cell 2 states 3,3,3,3
	m <- cbind(`epsilon[1]` = c(1, 1, 2, 3), `epsilon[2]` = c(3, 3, 3, 3),
	           `theta[1]` = rep(0.3, 4), `theta[2]` = rep(0.3, 4), `theta[3]` = rep(0.4, 4))
	out <- cell_prob(m, obj)
	expect_equal(unname(out), matrix(c(0.5, 0.25, 0.25, 0, 0, 1), 3, 2))
	expect_identical(colnames(out), c("epsilon[1]", "epsilon[2]"))
	expect_equal(unname(colSums(out)), c(1, 1))
})

test_that("cell_prob matches a reference implementation for i3 and i6 on random draws", {
	set.seed(11)
	for (K in c(3, 6)) {
		obj <- bn_manual_obj(list(), numeric(0), HMM_type = if (K == 3) "i3" else "i6")
		eps <- matrix(sample(seq_len(K), 50 * 5, replace = TRUE), 50, 5)
		colnames(eps) <- paste0("epsilon[", 1:5, "]")
		th <- matrix(1 / K, 50, K, dimnames = list(NULL, paste0("theta[", seq_len(K), "]")))
		m <- cbind(eps, th)
		out <- cell_prob(m, obj)
		expect_equal(dim(out), c(K, 5L))
		expect_equal(unname(out), bn_ref_cell_prob(m, K))
		expect_equal(unname(colSums(out)), rep(1, 5))
		expect_true(all(out >= 0 & out <= 1))
	}
})

test_that("cell_prob reports zero for states never drawn and handles one cell", {
	obj <- bn_manual_obj(list(), numeric(0), HMM_type = "i6")
	m <- cbind(`epsilon` = c(2, 2, 5, 5, 5), `theta[1]` = rep(1/6, 5))
	out <- cell_prob(m, obj)
	expect_equal(dim(out), c(6L, 1L))
	expect_equal(as.numeric(out), c(0, 0.4, 0, 0, 0.6, 0))
})

# ----------------------------------------------------------- getProbabilities

test_that("getProbabilities combines chains and splits thetas from epsilons", {
	obj <- bn_manual_obj(list(), numeric(0), HMM_type = "i3")
	chain <- function(e1, e2, th) {
		m <- cbind(`epsilon[1]` = e1, `epsilon[2]` = e2,
		           `theta[1]` = th[1], `theta[2]` = th[2], `theta[3]` = th[3])
		m
	}
	# region 1: two chains of 2 draws; region 2: two chains of 3 draws, one cell
	r1 <- list(chain(c(1, 1), c(3, 3), c(0.5, 0.2, 0.3)),
	           chain(c(2, 3), c(3, 1), c(0.1, 0.2, 0.7)))
	r2 <- list(cbind(`epsilon` = c(2, 2, 2), `theta[1]` = 0.2, `theta[2]` = 0.6, `theta[3]` = 0.2),
	           cbind(`epsilon` = c(2, 2, 1), `theta[1]` = 0.2, `theta[2]` = 0.5, `theta[3]` = 0.3))
	res <- getProbabilities(obj, list(r1, r2))
	expect_length(res@cnv_probabilities, 2)
	expect_length(res@cell_probabilities, 2)
	expect_equal(dim(res@cnv_probabilities[[1]]), c(4L, 3L))
	expect_equal(dim(res@cnv_probabilities[[2]]), c(6L, 3L))
	expect_equal(unname(res@cnv_probabilities[[1]][, 1]), c(0.5, 0.5, 0.1, 0.1))
	expect_equal(unname(colMeans(res@cnv_probabilities[[2]])), c(0.2, 0.55, 0.25))
	# cell 1: draws 1,1,2,3 -> (.5,.25,.25); cell 2: draws 3,3,3,1 -> (.25,0,.75)
	expect_equal(unname(res@cell_probabilities[[1]]), matrix(c(.5, .25, .25, .25, 0, .75), 3, 2))
	# region 2: draws 2,2,2,2,2,1 -> (1/6, 5/6, 0)
	expect_equal(as.numeric(res@cell_probabilities[[2]]), c(1/6, 5/6, 0))
	# the input object is not modified
	expect_length(obj@cnv_probabilities, 0)
})

# --------------------------------------------------------------- getGenesCells

test_that("getGenesCells maps regions to gene and cell indices", {
	m <- matrix(0, 5, 6, dimnames = list(paste0("g", 1:5), paste0("c", 1:6)))
	obj <- methods::new("MCMC_inferCNV")
	obj@expr.data <- m
	pred <- data.frame(cell_group_name = c("A", "A", "B", "B", "B"),
	                   gene_region_name = c("R1", "R1", "R2", "R2", "R2"),
	                   state = c(1, 1, 3, 3, 3),
	                   gene = c("g4", "g2", "g5", "g1", "g3"),
	                   stringsAsFactors = TRUE)
	groups <- data.frame(cell_group_name = c("A", "A", "B", "B", "C"),
	                     cell = c("c6", "c2", "c3", "c4", "c1"))
	obj@cnv_regions <- unique(pred$gene_region_name)
	res <- getGenesCells(obj, pred, groups)
	cg <- cellGene(res)
	expect_length(cg, 2)
	expect_equal(as.character(cg[[1]]$cnv_regions), "R1")
	expect_equal(cg[[1]]$Genes, c(2L, 4L))      # sorted by position, not by file order
	expect_equal(cg[[1]]$Cells, c(2L, 6L))
	expect_equal(cg[[1]]$State, 1)
	expect_equal(as.character(cg[[2]]$cnv_regions), "R2")
	expect_equal(cg[[2]]$Genes, c(1L, 3L, 5L))
	expect_equal(cg[[2]]$Cells, c(3L, 4L))
	expect_equal(cg[[2]]$State, 3)
})

# ------------------------------------------------------------ initializeObject

test_that("initializeObject builds regions, cell groups, args and model from files", {
	obj <- bn_make_obj(set_mean_sd = FALSE)
	expect_s4_class(obj, "MCMC_inferCNV")
	expect_equal(as.character(obj@cnv_regions), c("chr1-region_1", "chr1-region_2"))
	cg <- cellGene(obj)
	expect_length(cg, 2)
	expect_equal(cg[[1]]$Genes, 1:6)
	expect_equal(cg[[1]]$Cells, 5:6)
	expect_equal(cg[[1]]$State, 1)
	expect_equal(cg[[2]]$Genes, 7:12)
	expect_equal(cg[[2]]$Cells, 7:8)
	expect_equal(cg[[2]]$State, 3)
	expect_identical(obj@group_id, c(NA, NA, NA, NA, 1L, 1L, 2L, 2L))
	expect_identical(modelFile(obj), getArgs(obj)$model_file)
	expect_true(grepl("BUGS_Mixture_Model_i3$", modelFile(obj)))
	expect_equal(expirData(obj), bn_make_expr())
	expect_equal(obj@expr.data, bn_make_expr())
})

test_that("initializeObject honours resume_file_token in the file names", {
	dir <- tempfile("bn_token")
	bn_write_fixture_dir(dir)
	file.rename(file.path(dir, "17_HMM_pred.cell_groupings"), file.path(dir, "17_HMM_predTOK.cell_groupings"))
	file.rename(file.path(dir, "17_HMM_pred.pred_cnv_genes.dat"), file.path(dir, "17_HMM_predTOK.pred_cnv_genes.dat"))
	obj <- methods::new("MCMC_inferCNV")
	suppressMessages(capture.output(
		obj <- initializeObject(obj, bn_args(dir, tempfile("bn_out"), resume_file_token = "TOK"), bn_make_infercnv())))
	expect_length(cellGene(obj), 2)
	expect_equal(cellGene(obj)[[1]]$Cells, 5:6)
})

# ---------------------------------------------------------------------- MeanSD

test_that("MeanSD (i3) sets means mu-delta, mu, mu+delta and precision 1/sigma^2", {
	obj <- bn_make_obj(set_mean_sd = TRUE)
	ref <- bn_make_expr()[, 1:4]
	mu <- mean(ref)
	sigma <- sd(ref)
	delta <- abs(qnorm(0.05, 0, sigma))
	expect_equal(mu, 1)
	expect_equal(obj@mu, c(mu - delta, mu, mu + delta))
	expect_equal(obj@sig, rep(1 / sigma^2, 3))
})

test_that("MeanSD (i3) uses the observations when there are no references", {
	m <- bn_make_expr()
	inf <- bn_make_infercnv(m)
	obj <- bn_make_obj(m, set_mean_sd = FALSE)
	obj@reference_grouped_cell_indices <- list()
	suppressMessages(capture.output(obj <- MeanSD(obj, NULL, inf)))
	tum <- m[, 5:8]
	expect_equal(obj@mu[2], mean(tum))
	expect_equal(obj@sig[1], 1 / sd(tum)^2)
})

test_that("MeanSD (i6) takes per-level mean and precision from the hidden spike", {
	chr_info <- .get_hspike_chr_info(1, 1)
	chrs <- vapply(chr_info, function(x) x$name, "")
	cnvs <- vapply(chr_info, function(x) x$cnv, 0)
	genes_per_chr <- 3
	gene_chr <- rep(chrs, each = genes_per_chr)
	ng <- length(gene_chr)
	nc <- 4
	# deterministic values that differ between levels
	m <- matrix(0, ng, nc, dimnames = list(paste0("s", seq_len(ng)), paste0("sc", seq_len(nc))))
	for (i in seq_len(ng)) {
		lvl <- cnvs[match(gene_chr[i], chrs)]
		for (j in seq_len(nc)) m[i, j] <- lvl + 0.01 * i + 0.003 * j^2 * (-1)^(i + j)
	}
	go <- data.frame(chr = gene_chr, start = seq_len(ng), stop = seq_len(ng), row.names = rownames(m))
	hs <- methods::new("infercnv", expr.data = m, count.data = m, gene_order = go,
	                   reference_grouped_cell_indices = list(),
	                   observation_grouped_cell_indices = list(spike = seq_len(nc)))
	obj <- bn_sim_obj(numeric(0), numeric(0), "i6")
	obj@.hspike <- hs
	res <- MeanSD(obj, HMM_states = NULL, infercnv_obj = NULL)
	levels_sorted <- sort(unique(sprintf("cnv:%g", cnvs)))
	expect_length(res@mu, 6)
	expect_length(res@sig, 6)
	exp_mu <- exp_sig <- numeric(0)
	for (lv in levels_sorted) {
		v <- c(m[gene_chr %in% chrs[sprintf("cnv:%g", cnvs) == lv], ])
		exp_mu <- c(exp_mu, mean(v))
		exp_sig <- c(exp_sig, 1 / sd(v)^2)
	}
	expect_equal(unname(res@mu), exp_mu)
	expect_equal(unname(res@sig), exp_sig)
	expect_equal(levels_sorted, c("cnv:0.01", "cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3"))
	expect_true(all(diff(res@mu) > 0))
})

test_that(".get_gene_expr_by_cnv pools chromosomes with the same level and mean/sd helper is exact", {
	hs <- NULL
	chr_info <- .get_hspike_chr_info(1, 1)
	chrs <- vapply(chr_info, function(x) x$name, "")
	ng <- length(chrs)
	m <- matrix(seq_len(ng * 2), ng, 2, dimnames = list(paste0("s", seq_len(ng)), c("a", "b")))
	go <- data.frame(chr = chrs, start = seq_len(ng), stop = seq_len(ng), row.names = rownames(m))
	hs <- methods::new("infercnv", expr.data = m, count.data = m, gene_order = go,
	                   reference_grouped_cell_indices = list(),
	                   observation_grouped_cell_indices = list(spike = 1:2))
	by <- .get_gene_expr_by_cnv(hs)
	expect_setequal(names(by), c("cnv:0.01", "cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3"))
	# level 1 = chrA, chr_B, chr_C, chr_D, chr_E, chr_F : 6 genes x 2 cells
	idx1 <- which(chrs %in% c("chrA", "chr_B", "chr_C", "chr_D", "chr_E", "chr_F"))
	expect_equal(sort(unname(by[["cnv:1"]])), sort(c(m[idx1, ])))
	expect_equal(unname(by[["cnv:0.01"]]), unname(c(m[chrs == "chr_0", ])))
	ms <- .get_gene_expr_mean_sd_by_cnv(list(x = c(1, 2, 3, 4), y = c(5, 5)))
	expect_equal(ms$x$mean, 2.5)
	expect_equal(ms$x$sd, sqrt(5 / 3))
	expect_equal(ms$y$mean, 5)
	expect_equal(ms$y$sd, 0)
})

# ------------------------------------------------------------- run_gibb_sampling

bn_check_samples <- function(s, n_chains, n_cells, K) {
	expect_s3_class(s, "mcmc.list")
	expect_length(s, n_chains)
	for (ch in seq_len(n_chains)) {
		expect_equal(dim(s[[ch]]), c(1000L, n_cells + K))
		th <- s[[ch]][, grepl("theta", colnames(s[[ch]])), drop = FALSE]
		ep <- s[[ch]][, grepl("epsilon", colnames(s[[ch]])), drop = FALSE]
		expect_equal(ncol(th), K)
		expect_equal(ncol(ep), n_cells)
		expect_equal(unname(rowSums(th)), rep(1, 1000))
		expect_true(all(th >= 0 & th <= 1))
		expect_true(all(ep %in% seq_len(K)))
	}
}

test_that("run_gibb_sampling (i3) returns 3 chains of 1000 draws with valid theta and epsilon", {
	bn_skip_if_no_jags()
	obj <- bn_sim_obj(c(0.5, 1, 1.5), rep(400, 3), "i3")
	g <- matrix(c(0.5, 0.5, 0.5, 1, 1, 1, 1.5, 1.5, 1.5, 1, 1, 1), 3, 4)
	s <- run_gibb_sampling(g, obj)
	bn_check_samples(s, 3, 4, 3)
	expect_identical(colnames(s[[1]]),
	                 c(paste0("epsilon[", 1:4, "]"), paste0("theta[", 1:3, "]")))
})

test_that("run_gibb_sampling (i3) assigns separated cells to the state with the nearest mean", {
	bn_skip_if_no_jags()
	obj <- bn_sim_obj(c(0.5, 1, 1.5), rep(400, 3), "i3")
	# cells: state 1, state 1, state 1, state 3, state 3 (8 genes each, exact means)
	g <- matrix(rep(c(0.5, 0.5, 0.5, 1.5, 1.5), each = 8), 8, 5)
	s <- run_gibb_sampling(g, obj)
	bn_check_samples(s, 3, 5, 3)
	cp <- cell_prob(do.call(rbind, s), obj)
	expect_equal(dim(cp), c(3L, 5L))
	expect_true(all(cp[1, 1:3] > 0.99))
	expect_true(all(cp[3, 4:5] > 0.99))
	# conjugate check: theta | eps ~ Dirichlet(1 + counts) = Dirichlet(4, 1, 3),
	# posterior mean (0.5, 0.125, 0.375); 3000 draws, tolerance 0.05
	th <- colMeans(cnv_prob(do.call(rbind, s)))
	expect_equal(unname(th), c(0.5, 0.125, 0.375), tolerance = 0.1)
	expect_true(all(abs(th - c(0.5, 0.125, 0.375)) < 0.05))
})

test_that("run_gibb_sampling (i6) returns 6 chains and recovers states 2 and 5", {
	bn_skip_if_no_jags()
	obj <- bn_sim_obj(1:6, rep(400, 6), "i6")
	g <- matrix(rep(c(2, 2, 5), each = 5), 5, 3)
	s <- run_gibb_sampling(g, obj)
	bn_check_samples(s, 6, 3, 6)
	cp <- cell_prob(do.call(rbind, s), obj)
	expect_equal(dim(cp), c(6L, 3L))
	expect_true(all(cp[2, 1:2] > 0.99))
	expect_true(cp[5, 3] > 0.99)
	# Dirichlet(1 + counts) = (1,3,1,1,2,1)/9 posterior mean
	th <- colMeans(cnv_prob(do.call(rbind, s)))
	expect_true(all(abs(th - c(1, 3, 1, 1, 2, 1) / 9) < 0.05))
})

test_that("run_gibb_sampling accepts a data frame and a single-cell vector", {
	bn_skip_if_no_jags()
	obj <- bn_sim_obj(c(0.5, 1, 1.5), rep(400, 3), "i3")
	df <- as.data.frame(matrix(rep(c(0.5, 1.5), each = 6), 6, 2))
	s <- run_gibb_sampling(df, obj)
	bn_check_samples(s, 3, 2, 3)
	cp <- cell_prob(do.call(rbind, s), obj)
	expect_true(cp[1, 1] > 0.99 && cp[3, 2] > 0.99)
	# one cell, six genes: vector input is wrapped into a one-column data frame
	s1 <- run_gibb_sampling(rep(1.5, 6), obj)
	expect_length(s1, 3)
	expect_equal(dim(s1[[1]]), c(1000L, 4L))
	expect_true(all(s1[[1]][, "epsilon"] == 3))
})

test_that("posterior concentrates on state 2 for data at the middle mean regardless of precision scale", {
	bn_skip_if_no_jags()
	obj <- bn_sim_obj(c(0.5, 1, 1.5), rep(100, 3), "i3")
	g <- matrix(1, 10, 3)
	s <- run_gibb_sampling(g, obj)
	cp <- cell_prob(do.call(rbind, s), obj)
	expect_true(all(cp[2, ] > 0.99))
	th <- colMeans(cnv_prob(do.call(rbind, s)))
	# Dirichlet(1,4,1) mean (1/6, 4/6, 1/6)
	expect_true(all(abs(th - c(1, 4, 1) / 6) < 0.05))
})

# ------------------------------------------------- nonParallel / withParallel / runMCMC

test_that("nonParallel returns one mcmc.list per region and list(NULL) for an empty one", {
	bn_skip_if_no_jags()
	obj <- bn_make_obj()
	obj@cell_gene[[3]] <- list(cnv_regions = factor("empty"), Genes = 1:2, Cells = integer(0), State = 2)
	res <- nonParallel(obj)
	expect_length(res, 3)
	expect_s3_class(res[[1]], "mcmc.list")
	expect_equal(dim(res[[1]][[1]]), c(1000L, 2L + 3L))
	expect_equal(dim(res[[2]][[1]]), c(1000L, 2L + 3L))
	expect_identical(res[[3]], list(NULL))
})

test_that("withParallel with 2 cores gives results shaped like nonParallel", {
	bn_skip_if_no_jags()
	testthat::skip_on_os("windows")
	obj <- bn_make_obj()
	obj@args$CORES <- 2
	res <- withParallel(obj)
	expect_length(res, 2)
	for (r in res) {
		expect_s3_class(r, "mcmc.list")
		expect_length(r, 3)
		expect_equal(dim(r[[1]]), c(1000L, 5L))
	}
	cp <- cell_prob(do.call(rbind, res[[1]]), obj)
	expect_true(all(cp[3, ] > 0.99))
})

test_that("runMCMC fills cnv and cell probabilities that identify the true states", {
	bn_skip_if_no_jags()
	obj <- bn_make_obj()
	expect_equal(unname(obj@mu[2]), 1)
	res <- runMCMC(obj, FALSE)
	expect_length(res@cnv_probabilities, 2)
	expect_length(res@cell_probabilities, 2)
	for (i in 1:2) {
		expect_equal(dim(res@cnv_probabilities[[i]]), c(3000L, 3L))
		expect_equal(unname(rowSums(res@cnv_probabilities[[i]])), rep(1, 3000))
		expect_equal(dim(res@cell_probabilities[[i]]), c(3L, 2L))
		expect_equal(unname(colSums(res@cell_probabilities[[i]])), c(1, 1))
	}
	# region 1 (amplified data) is state 3, region 2 (reference level) is state 2
	expect_true(all(res@cell_probabilities[[1]][3, ] > 0.99))
	expect_true(all(res@cell_probabilities[[2]][2, ] > 0.99))
	# two cells in one state, Dirichlet(1,1,1) prior: posterior mean (1,1,3)/5 ordering
	m1 <- colMeans(res@cnv_probabilities[[1]])
	m2 <- colMeans(res@cnv_probabilities[[2]])
	expect_true(abs(m1[3] - 0.6) < 0.06 && abs(m1[1] - 0.2) < 0.06 && abs(m1[2] - 0.2) < 0.06)
	expect_true(abs(m2[2] - 0.6) < 0.06 && abs(m2[1] - 0.2) < 0.06 && abs(m2[3] - 0.2) < 0.06)
	# the sampler's other slots are kept
	expect_equal(res@mu, obj@mu)
	expect_equal(cellGene(res), cellGene(obj))
})

# ----------------------------------------------------------------- reassignCNV

test_that("reassignCNV moves each CNV to its highest posterior state and rewrites HMM_states", {
	# means: r1 (.1,.2,.7) HMM 1 -> 3 ; r2 (.2,.6,.2) HMM 2 -> 2 ; r3 (.5,.3,.2) HMM 3 -> 1
	obj <- bn_manual_obj(list(bn_theta_draws(c(.1, .2, .7)),
	                          bn_theta_draws(c(.2, .6, .2)),
	                          bn_theta_draws(c(.5, .3, .2))),
	                     state = c(1, 2, 3))
	H <- matrix(0, 6, 6)
	H[1:2, 1:2] <- 1; H[3:4, 3:4] <- 2; H[5:6, 5:6] <- 3
	out <- suppressMessages(reassignCNV(obj, H))
	expect_length(out, 2)
	expect_equal(vapply(cellGene(out[[1]]), function(x) x$State, 0), c(3, 2, 1))
	expectedH <- matrix(0, 6, 6)
	expectedH[1:2, 1:2] <- 3; expectedH[3:4, 3:4] <- 2; expectedH[5:6, 5:6] <- 1
	expect_equal(out[[2]], expectedH)
	# nothing else on the object changed
	expect_equal(out[[1]]@cnv_probabilities, obj@cnv_probabilities)
	expect_equal(lapply(cellGene(out[[1]]), function(x) x[c("cnv_regions", "Genes", "Cells")]),
	             lapply(cellGene(obj), function(x) x[c("cnv_regions", "Genes", "Cells")]))
})

test_that("reassignCNV with 6 states picks the argmax over all six", {
	p <- function(v) { m <- matrix(rep(v, each = 4), 4); colnames(m) <- paste0("theta[", 1:6, "]"); m }
	obj <- bn_manual_obj(list(p(c(.05, .05, .1, .1, .6, .1)), p(c(.5, .1, .1, .1, .1, .1 + 0.0))),
	                     state = c(3, 3), HMM_type = "i6")
	# second vector sums to 1.0: (.5, .1 x5)
	H <- matrix(0, 4, 4)
	out <- suppressMessages(reassignCNV(obj, H))
	expect_equal(vapply(cellGene(out[[1]]), function(x) x$State, 0), c(5, 1))
	expect_equal(out[[2]][1:2, 1:2], matrix(5, 2, 2))
	expect_equal(out[[2]][3:4, 3:4], matrix(1, 2, 2))
	expect_equal(sum(out[[2]]), 5 * 4 + 1 * 4)
})

# ------------------------------------------------------------------- removeCNV

bn_remove_obj <- function() {
	# region means of P(normal = state 2): .2, .6, .1 ; threshold 0.5 removes region 2
	bn_manual_obj(list(bn_theta_draws(c(.5, .2, .3)),
	                   bn_theta_draws(c(.2, .6, .2)),
	                   bn_theta_draws(c(.1, .1, .8))),
	              state = c(1, 3, 3))
}

test_that("removeCNV drops CNVs above the normal-probability threshold and writes the table", {
	obj <- bn_remove_obj()
	H <- matrix(0, 6, 6)
	H[1:2, 1:2] <- 1; H[3:4, 3:4] <- 3; H[5:6, 5:6] <- 3
	out <- suppressMessages(removeCNV(obj, H))
	o <- out[[1]]
	expect_length(o@cell_gene, 2)
	expect_equal(vapply(o@cell_gene, function(x) as.character(x$cnv_regions), ""), c("r1", "r3"))
	expect_length(o@cnv_probabilities, 2)
	expect_length(o@cell_probabilities, 2)
	expect_equal(o@cell_probabilities, list(1L, 3L))
	expect_equal(unname(colMeans(o@cnv_probabilities[[2]])), c(.1, .1, .8))
	# removed CNV's block is set to the normal state (2); the rest untouched
	expectedH <- H
	expectedH[3:4, 3:4] <- 2
	expect_equal(out[[2]], expectedH)
	# table of the remaining CNVs
	f <- file.path(getArgs(o)$out_dir, "CNV_State_Probabilities.dat")
	expect_true(file.exists(f))
	tab <- utils::read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
	expect_equal(colnames(tab), c("r1", "r3"))
	expect_equal(rownames(tab), c("State:1", "State:2", "State:3"))
	expect_equal(as.matrix(tab), matrix(c(.5, .2, .3, .1, .1, .8), 3, 2, dimnames = list(c("State:1", "State:2", "State:3"), c("r1", "r3"))),
	             tolerance = 1e-8)
})

test_that("removeCNV keeps everything and writes nothing when no CNV exceeds the threshold", {
	obj <- bn_remove_obj()
	obj@args$BayesMaxPNormal <- 0.9
	H <- matrix(1, 6, 6)
	out <- suppressMessages(removeCNV(obj, H))
	expect_identical(out[[1]], obj)
	expect_identical(out[[2]], H)
	expect_false(dir.exists(getArgs(obj)$out_dir))
})

test_that("removeCNV removes every CNV above a strict threshold using >, not >=", {
	obj <- bn_remove_obj()
	# region 2 has P(normal) = 0.6 up to floating point: threshold 0.6 must keep it if the mean is <= 0.6
	thr <- colMeans(obj@cnv_probabilities[[2]])[2]
	obj@args$BayesMaxPNormal <- thr
	out <- suppressMessages(removeCNV(obj, matrix(0, 6, 6)))
	expect_length(out[[1]]@cell_gene, 3)
	obj@args$BayesMaxPNormal <- thr - 1e-6
	out <- suppressMessages(removeCNV(obj, matrix(0, 6, 6)))
	expect_length(out[[1]]@cell_gene, 2)
})

test_that("removeCNV uses state 3 as normal and State:1..6 row names for i6", {
	p <- function(v) { m <- matrix(rep(v, each = 4), 4); colnames(m) <- paste0("theta[", 1:6, "]"); m }
	obj <- bn_manual_obj(list(p(c(.1, .1, .6, .1, .05, .05)),
	                          p(c(.1, .1, .1, .1, .5, .1)),
	                          p(c(.3, .3, .1, .1, .1, .1))),
	                     state = c(5, 5, 1), HMM_type = "i6")
	out <- suppressMessages(removeCNV(obj, matrix(0, 6, 6)))
	expect_equal(vapply(out[[1]]@cell_gene, function(x) as.character(x$cnv_regions), ""), c("r2", "r3"))
	expect_equal(out[[2]][1:2, 1:2], matrix(3, 2, 2))
	tab <- utils::read.table(file.path(getArgs(out[[1]])$out_dir, "CNV_State_Probabilities.dat"),
	                         header = TRUE, sep = "\t", check.names = FALSE)
	expect_equal(rownames(tab), paste0("State:", 1:6))
	expect_equal(colnames(tab), c("r2", "r3"))
})

# ----------------------------------------------------------------- removeCells

test_that("removeCells returns the object unchanged when no cell exceeds the threshold", {
	cp <- list(matrix(c(.8, .1, .1, .7, .2, .1), 3, 2), matrix(c(.1, .3, .6, .2, .3, .5), 3, 2))
	obj <- bn_manual_obj(list(bn_theta_draws(c(.4, .3, .3)), bn_theta_draws(c(.4, .3, .3))),
	                     state = c(1, 3), cell_probs = cp)
	out <- removeCells(obj, matrix(0, 4, 4))
	expect_identical(out, obj)
})

# ----------------------------------------------------------- filterHighPNormals

test_that("filterHighPNormals with no post-MCMC method only records the threshold", {
	obj <- bn_remove_obj()
	obj@args$postMcmcMethod <- NULL
	H <- matrix(2, 6, 6)
	res <- filterHighPNormals(obj, H, BayesMaxPNormal = 0.3, useRaster = FALSE)
	expect_length(res, 2)
	expect_equal(getArgs(res[[1]])$BayesMaxPNormal, 0.3)
	expect_identical(res[[2]], H)
	expect_length(cellGene(res[[1]]), 3)
})

test_that("filterHighPNormals removes high-normal CNVs then reassigns the rest", {
	obj <- bn_manual_obj(list(bn_theta_draws(c(.1, .2, .7)),      # HMM 1 -> reassigned to 3
	                          bn_theta_draws(c(.2, .6, .2)),      # removed, set to normal 2
	                          bn_theta_draws(c(.6, .2, .2))),     # HMM 3 -> reassigned to 1
	                     state = c(1, 3, 3), BayesMaxPNormal = 0)
	H <- matrix(0, 6, 6)
	H[1:2, 1:2] <- 1; H[3:4, 3:4] <- 3; H[5:6, 5:6] <- 3
	res <- suppressMessages(filterHighPNormals(obj, H, BayesMaxPNormal = 0.5, useRaster = FALSE))
	o <- res[[1]]
	expect_equal(getArgs(o)$BayesMaxPNormal, 0.5)
	expect_equal(vapply(cellGene(o), function(x) as.character(x$cnv_regions), ""), c("r1", "r3"))
	expect_equal(vapply(cellGene(o), function(x) x$State, 0), c(3, 1))
	expectedH <- matrix(0, 6, 6)
	expectedH[1:2, 1:2] <- 3; expectedH[3:4, 3:4] <- 2; expectedH[5:6, 5:6] <- 1
	expect_equal(res[[2]], expectedH)
})

test_that("filterHighPNormals without reassignment keeps the HMM states of kept CNVs", {
	obj <- bn_manual_obj(list(bn_theta_draws(c(.1, .2, .7)),
	                          bn_theta_draws(c(.2, .6, .2)),
	                          bn_theta_draws(c(.6, .2, .2))),
	                     state = c(1, 3, 3), BayesMaxPNormal = 0, reassignCNVs = FALSE)
	H <- matrix(0, 6, 6)
	H[1:2, 1:2] <- 1; H[3:4, 3:4] <- 3; H[5:6, 5:6] <- 3
	res <- suppressMessages(filterHighPNormals(obj, H, BayesMaxPNormal = 0.5, useRaster = FALSE))
	expect_equal(vapply(cellGene(res[[1]]), function(x) x$State, 0), c(1, 3))
	expectedH <- H
	expectedH[3:4, 3:4] <- 2
	expect_equal(res[[2]], expectedH)
})

# ------------------------------------------------------------ inferCNVBayesNet

test_that("inferCNVBayesNet rejects a missing model file and too many cores", {
	dir <- tempfile("bn_err")
	bn_write_fixture_dir(dir)
	expect_error(
		suppressMessages(inferCNVBayesNet(file_dir = dir, infercnv_obj = bn_make_infercnv(),
		                                  HMM_states = NULL, out_dir = tempfile("bn_o"),
		                                  resume_file_token = "",
		                                  model_file = file.path(tempdir(), "does_not_exist"),
		                                  HMM_type = "i3", k_obs_groups = 1, cluster_by_groups = FALSE,
		                                  no_plot = TRUE, useRaster = FALSE)),
		"Cannot find the model file")
	testthat::skip_if(is.na(parallel::detectCores()))
	expect_error(
		suppressMessages(inferCNVBayesNet(file_dir = dir, infercnv_obj = bn_make_infercnv(),
		                                  HMM_states = NULL, out_dir = tempfile("bn_o"),
		                                  resume_file_token = "", CORES = parallel::detectCores() + 1000,
		                                  HMM_type = "i3", k_obs_groups = 1, cluster_by_groups = FALSE,
		                                  no_plot = TRUE, useRaster = FALSE)),
		"Too many cores")
})

test_that("inferCNVBayesNet (i3) runs end to end, saves the object and applies no_plot", {
	bn_skip_if_no_jags()
	dir <- tempfile("bn_e2e")
	bn_write_fixture_dir(dir)
	out_dir <- file.path(tempfile("bn_e2e_out"))
	inf <- bn_make_infercnv()
	res <- suppressMessages(capture.output(
		obj <- inferCNVBayesNet(file_dir = dir, infercnv_obj = inf, HMM_states = NULL,
		                        out_dir = out_dir, resume_file_token = "",
		                        postMcmcMethod = "removeCNV", plotingProbs = TRUE,
		                        diagnostics = TRUE, HMM_type = "i3", k_obs_groups = 1,
		                        cluster_by_groups = FALSE, no_plot = TRUE, useRaster = FALSE)))
	expect_s4_class(obj, "MCMC_inferCNV")
	expect_true(dir.exists(out_dir))
	expect_true(file.exists(file.path(out_dir, "MCMC_inferCNV_obj.rds")))
	saved <- readRDS(file.path(out_dir, "MCMC_inferCNV_obj.rds"))
	expect_equal(saved, obj)
	# no_plot overrides plotingProbs / diagnostics and no plot files appear
	a <- getArgs(obj)
	expect_false(a$plotingProbs)
	expect_false(a$diagnostics)
	expect_equal(a$BayesMaxPNormal, 0)
	expect_equal(a$CORES, 1)
	expect_true(grepl("BUGS_Mixture_Model_i3$", a$model_file))
	expect_length(list.files(out_dir, pattern = "pdf$"), 0)
	# state assignment and probabilities
	sigma <- sd(bn_make_expr()[, 1:4])
	expect_equal(obj@mu, c(1 - abs(qnorm(0.05, 0, sigma)), 1, 1 + abs(qnorm(0.05, 0, sigma))))
	expect_equal(obj@sig, rep(1 / sigma^2, 3))
	expect_length(obj@cnv_probabilities, 2)
	expect_true(all(obj@cell_probabilities[[1]][3, ] > 0.99))
	expect_true(all(obj@cell_probabilities[[2]][2, ] > 0.99))
})

# -------------------------------------------------------------------- plotting

test_that("plot_cell_prob and plot_cnv_prob build ggplots from the probability tables", {
	testthat::skip_if_not_installed("ggplot2")
	cp <- matrix(c(.5, .25, .25, 0, 0, 1), 3, 2, dimnames = list(NULL, c("epsilon[1]", "epsilon[2]")))
	p <- suppressWarnings(plot_cell_prob(as.data.frame(cp), "region X", "i3"))
	expect_s3_class(p, "ggplot")
	expect_equal(p$labels$title, "region X")
	d <- p$data
	expect_equal(nrow(d), 6)
	expect_equal(as.character(d$mag), rep(as.character(1:3), 2))
	expect_equal(d$value, c(.5, .25, .25, 0, 0, 1))
	p6 <- suppressWarnings(plot_cell_prob(as.data.frame(matrix(1 / 6, 6, 2)), "t", "i6"))
	expect_equal(nlevels(p6$data$mag), 6)
	th <- bn_theta_draws(c(.2, .3, .5))
	pc <- suppressWarnings(plot_cnv_prob(th, "cnv Y", "i3"))
	expect_s3_class(pc, "ggplot")
	expect_equal(nrow(pc$data), 12)
	expect_equal(pc$labels$title, "cnv Y")
	expect_equal(sort(unique(as.character(pc$data$State))), c("1", "2", "3"))
	expect_equal(sort(pc$data$Probability), sort(c(th)))
})

test_that("plotProbabilities writes the cell and CNV PDFs named by the threshold", {
	testthat::skip_if_not_installed("ggplot2")
	cp <- matrix(c(.5, .25, .25, 0, 0, 1), 3, 2)
	obj <- bn_manual_obj(list(bn_theta_draws(c(.4, .3, .3))), state = 1,
	                     cell_probs = list(cp), plotingProbs = TRUE)
	dir.create(getArgs(obj)$out_dir)
	suppressWarnings(suppressMessages(plotProbabilities(obj)))
	expect_true(file.exists(file.path(getArgs(obj)$out_dir, "cellProbs.0.5.pdf")))
	expect_true(file.exists(file.path(getArgs(obj)$out_dir, "cnvProbs.0.5.pdf")))
	# nothing is written when plotingProbs is FALSE
	obj2 <- bn_manual_obj(list(bn_theta_draws(c(.4, .3, .3))), state = 1, cell_probs = list(cp))
	dir.create(getArgs(obj2)$out_dir)
	suppressMessages(plotProbabilities(obj2))
	expect_length(list.files(getArgs(obj2)$out_dir), 0)
})
