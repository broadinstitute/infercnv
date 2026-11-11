# Unit tests for the HMM state calling (R/inferCNV_HMM.R, R/inferCNV_i3HMM.R)
# and the CNV region reporting that consumes the HMM states.
#
# All fixtures are tiny synthetic matrices with a block of genes that is
# clearly lost or gained; expected states and regions are derived by hand.
# Helper names carry the prefix hr_ so they cannot clash with sibling files.

hr_old_threshold <- futile.logger::flog.threshold()
futile.logger::flog.threshold(futile.logger::FATAL)
withr::defer(futile.logger::flog.threshold(hr_old_threshold), testthat::teardown_env())

## ---------------------------------------------------------------------------
## fixtures
## ---------------------------------------------------------------------------

## six-state (i6) emission parameters: state i has mean hr_means6[i], shared sd
hr_means6 <- c(0, 0.5, 1, 1.5, 2, 3)

hr_cnv_mean_sd <- function(sd=0.1) {
	lv <- c("cnv:0.01", "cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3")
	out <- lapply(seq_along(lv), function(i) list(mean=hr_means6[i], sd=sd))
	names(out) <- lv
	out
}

## sd(n cells) = 0.2 / sqrt(n) exactly, as a list of lm fits keyed by cnv level,
## in the form the package builds (response log(sd), predictor log(num_cells)).
hr_sd_fits <- function() {
	d <- data.frame(num_cells=1:10)
	d$y <- log(0.2) - 0.5 * log(d$num_cells)
	fit <- lm(y ~ log(num_cells), data=d)
	lv <- c("cnv:0.01", "cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3")
	out <- rep(list(fit), length(lv))
	names(out) <- lv
	out
}

## chr c1 has 36 genes, chr c2 has 12 genes; 12 cells.
## cells 1-4 reference ("ref"), 5-8 tumor group A, 9-12 tumor group B.
## Group A has a loss (value 0.5) over genes 13-24 of c1; all other values are 1.
## Neutral flanks are 12 genes long so the HMM switches state cleanly.
hr_make_obj <- function(subclusters=NULL) {
	genes <- paste0("g", 1:48)
	ge <- data.frame(chr=rep(c("c1", "c2"), c(36, 12)), start=1:48, stop=1:48 + 1L,
					 stringsAsFactors=FALSE)
	rownames(ge) <- genes
	m <- matrix(1, nrow=48, ncol=12, dimnames=list(genes, paste0("cell", 1:12)))
	m[13:24, 5:8] <- 0.5
	new("infercnv", expr.data=m, count.data=m, gene_order=ge,
		reference_grouped_cell_indices=list(ref=1:4),
		observation_grouped_cell_indices=list(A=5:8, B=9:12),
		tumor_subclusters=subclusters, options=list(), .hspike=NULL)
}

## an object whose expr.data already holds HMM states, for the region functions
hr_make_state_obj <- function(subclusters=NULL) {
	ge <- data.frame(chr=c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr2", "chr3"),
					 start=c(100, 200, 300, 400, 500, 10, 20, 30, 5),
					 stop=c(150, 250, 350, 450, 550, 15, 25, 35, 9),
					 stringsAsFactors=FALSE)
	rownames(ge) <- paste0("gene", 1:9)
	## columns: ref1, ref2, tumA1, tumA2, tumA3
	m <- cbind(ref1=c(3, 3, 3, 3, 3, 3, 3, 3, 3),
			   ref2=c(3, 3, 3, 3, 3, 3, 3, 3, 3),
			   tumA1=c(3, 2, 2, 3, 3, 5, 5, 3, 1),
			   tumA2=c(3, 2, 2, 3, 3, 5, 5, 3, 1),
			   tumA3=c(3, 3, 3, 3, 4, 3, 3, 3, 1))
	rownames(m) <- rownames(ge)
	new("infercnv", expr.data=m, count.data=m, gene_order=ge,
		reference_grouped_cell_indices=list(ref=1:2),
		observation_grouped_cell_indices=list(tumA=3:5),
		tumor_subclusters=subclusters, options=list(), .hspike=NULL)
}

## hidden-spike object with the 11 fake chromosomes of .get_hspike_chr_info.
## Every chromosome has 2 genes and the 3 cells carry cnv_level * (1 + e) with
## e = (-.1, 0, .1) per gene, so the per-level sd is known in closed form.
hr_make_spike <- function() {
	info <- .get_hspike_chr_info(1, 1)
	chr <- unlist(lapply(info, function(i) rep(i$name, 2)))
	lvl <- unlist(lapply(info, function(i) rep(i$cnv, 2)))
	ge <- data.frame(chr=chr, start=seq_along(chr), stop=seq_along(chr), stringsAsFactors=FALSE)
	rownames(ge) <- paste0("s", seq_along(chr))
	m <- matrix(0, nrow=length(chr), ncol=3)
	e <- matrix(c(-1, 0, 1, -1, 0, 1), nrow=2, byrow=TRUE) * 0.01
	for (i in seq_along(info)) {
		m[(2 * i - 1):(2 * i), ] <- info[[i]]$cnv + e
	}
	rownames(m) <- rownames(ge)
	colnames(m) <- paste0("sc", 1:3)
	new("infercnv", expr.data=m, count.data=m, gene_order=ge,
		reference_grouped_cell_indices=list(),
		observation_grouped_cell_indices=list(spike=1:3),
		tumor_subclusters=NULL, options=list(), .hspike=NULL)
}

## ---------------------------------------------------------------------------
## spike distributions
## ---------------------------------------------------------------------------

test_that(".get_gene_expr_by_cnv pools values by cnv level in chromosome order", {
	spike <- hr_make_spike()
	by_cnv <- .get_gene_expr_by_cnv(spike)

	## order of first appearance in .get_hspike_chr_info
	expect_equal(names(by_cnv), c("cnv:1", "cnv:0.01", "cnv:0.5", "cnv:1.5", "cnv:2", "cnv:3"))

	## 6 neutral chromosomes (2 genes x 3 cells) pooled, 6 values per other level
	expect_length(by_cnv[["cnv:1"]], 36)
	expect_length(by_cnv[["cnv:0.5"]], 6)

	## values of a 2-gene chromosome are read column by column:
	## cell1 genes 1-2, cell2 genes 1-2, cell3 genes 1-2
	expect_equal(by_cnv[["cnv:0.5"]], 0.5 + c(-0.01, -0.01, 0, 0, 0.01, 0.01))
})

test_that(".get_gene_expr_mean_sd_by_cnv computes mean and sample sd per level", {
	by_cnv <- list("cnv:1"=c(1, 2, 3, 4), "cnv:3"=c(2, 4))
	out <- .get_gene_expr_mean_sd_by_cnv(by_cnv)
	expect_equal(names(out), c("cnv:1", "cnv:3"))
	## mean 2.5, sum of squared deviations 5, sample sd sqrt(5/3)
	expect_equal(out[["cnv:1"]]$mean, 2.5)
	expect_equal(out[["cnv:1"]]$sd, sqrt(5 / 3))
	## mean 3, sum sq 2, sd sqrt(2)
	expect_equal(out[["cnv:3"]]$mean, 3)
	expect_equal(out[["cnv:3"]]$sd, sqrt(2))
})

test_that("get_spike_dists gives the closed-form mean and sd for every level", {
	spike <- hr_make_spike()
	d <- get_spike_dists(spike)
	expect_equal(names(d), c("cnv:1", "cnv:0.01", "cnv:0.5", "cnv:1.5", "cnv:2", "cnv:3"))
	for (lv in c("cnv:0.01", "cnv:0.5", "cnv:1.5", "cnv:2", "cnv:3")) {
		lvl <- as.numeric(sub("cnv:", "", lv))
		expect_equal(d[[lv]]$mean, lvl, tolerance=1e-12)
		## deviations (-.01,-.01,0,0,.01,.01): sum of squares 4e-4, n-1 = 5
		expect_equal(d[[lv]]$sd, sqrt(4e-4 / 5), tolerance=1e-12)
	}
	## the neutral pool is 6 chromosomes of the same pattern: sum sq 6 * 4e-4, n-1 = 35
	expect_equal(d[["cnv:1"]]$mean, 1, tolerance=1e-12)
	expect_equal(d[["cnv:1"]]$sd, sqrt(6 * 4e-4 / 35), tolerance=1e-12)
})

test_that("get_spike_dists stops on a NULL spike object", {
	expect_error(get_spike_dists(NULL))
})

test_that("get_hspike_cnv_mean_sd_trend_by_num_cells_fit returns one lm per cnv level", {
	## the fitted slope is not asserted: it averages over the wrong margin (see the bug report)
	set.seed(42)
	fits <- get_hspike_cnv_mean_sd_trend_by_num_cells_fit(hr_make_spike())
	expect_equal(names(fits), c("cnv:1", "cnv:0.01", "cnv:0.5", "cnv:1.5", "cnv:2", "cnv:3"))
	for (lv in names(fits)) {
		expect_s3_class(fits[[lv]], "lm")
		expect_equal(names(coef(fits[[lv]])), c("(Intercept)", "log(num_cells)"))
		expect_true(all(is.finite(coef(fits[[lv]]))))
	}
})

## ---------------------------------------------------------------------------
## .get_HMM and emission parameters
## ---------------------------------------------------------------------------

test_that(".get_HMM builds the 6-state transition matrix, start vector and emissions", {
	t <- 0.01
	h <- .get_HMM(hr_cnv_mean_sd(0.2), t)
	expected_pi <- matrix(t, 6, 6)
	diag(expected_pi) <- 1 - 5 * t
	expect_equal(h$state_transitions, expected_pi)
	expect_equal(rowSums(h$state_transitions), rep(1, 6))
	expect_equal(h$delta, c(t, t, 1 - 5 * t, t, t, t))
	expect_equal(sum(h$delta), 1)
	expect_equal(h$state_emission_params$mean, hr_means6)
	expect_equal(h$state_emission_params$sd, rep(0.2, 6))
})

test_that(".get_HMM orders emission params by cnv level, not by list order", {
	msd <- hr_cnv_mean_sd(0.1)
	msd <- msd[c("cnv:3", "cnv:1", "cnv:0.01", "cnv:2", "cnv:0.5", "cnv:1.5")]
	msd[["cnv:2"]]$sd <- 0.7
	h <- .get_HMM(msd, 1e-6)
	expect_equal(h$state_emission_params$mean, hr_means6)
	expect_equal(h$state_emission_params$sd, c(0.1, 0.1, 0.1, 0.1, 0.7, 0.1))
})

test_that(".get_state_emission_params replaces sd with the fitted exp(a + b log n)", {
	fits <- hr_sd_fits()
	msd <- hr_cnv_mean_sd(9)  # sd values must be overwritten
	## sd(n) = 0.2 / sqrt(n): n=4 -> 0.1, n=1 -> 0.2, n=9 -> 0.2/3
	p4 <- .get_state_emission_params(4, msd, fits)
	expect_equal(p4$mean, hr_means6)
	expect_equal(p4$sd, rep(0.1, 6), tolerance=1e-9)
	expect_equal(.get_state_emission_params(1, msd, fits)$sd, rep(0.2, 6), tolerance=1e-9)
	expect_equal(.get_state_emission_params(9, msd, fits)$sd, rep(0.2 / 3, 6), tolerance=1e-9)
})

test_that(".get_state_emission_params uses a separate fit for each level", {
	fits <- hr_sd_fits()
	d <- data.frame(num_cells=1:10)
	d$y <- log(1.0) - 1 * log(d$num_cells)  # sd(n) = 1/n
	fits[["cnv:2"]] <- lm(y ~ log(num_cells), data=d)
	p <- .get_state_emission_params(4, hr_cnv_mean_sd(), fits)
	expect_equal(p$sd, c(0.1, 0.1, 0.1, 0.1, 0.25, 0.1), tolerance=1e-9)
})

test_that(".compare_obs_vs_fit_sc_cnv_sd compares observed sd to the fit at 1 cell", {
	fits <- hr_sd_fits()
	msd <- hr_cnv_mean_sd(0.3)
	df <- .compare_obs_vs_fit_sc_cnv_sd(NULL, cnv_mean_sd=msd, cnv_level_to_mean_sd_fit=fits)
	expect_equal(df$cnv, names(msd))
	expect_equal(df$sd_obs, rep(0.3, 6))
	expect_equal(unname(df$sd_pred), rep(0.2, 6), tolerance=1e-9)
})

## ---------------------------------------------------------------------------
## Viterbi
## ---------------------------------------------------------------------------

hr_dthmm <- function(x, sd=0.1, t=1e-6) {
	h <- .get_HMM(hr_cnv_mean_sd(sd), t)
	HiddenMarkov::dthmm(x=x, Pi=h$state_transitions, delta=h$delta, distn="norm",
						pm=h$state_emission_params)
}

## exhaustive search over all state paths using the scoring documented in
## Viterbi.dthmm.adj: per-gene emission weights e_s proportional to
## 1/(-log P(|Z| > |x-mu_s|/sd)), normalized to sum 1, plus log transition terms.
hr_bruteforce_viterbi <- function(x, means, sd, Pi, delta) {
	n <- length(x)
	m <- length(means)
	em <- t(sapply(x, function(v) {
		e <- 1 / (-pnorm(abs(v - means) / sd, log.p=TRUE, lower.tail=FALSE))
		log(e / sum(e))
	}))
	paths <- as.matrix(expand.grid(rep(list(seq_len(m)), n)))
	best <- -Inf
	best_path <- NULL
	for (r in seq_len(nrow(paths))) {
		p <- paths[r, ]
		s <- log(delta[p[1]]) + em[1, p[1]]
		for (i in seq_len(n)[-1]) {
			s <- s + log(Pi[p[i - 1], p[i]]) + em[i, p[i]]
		}
		if (s > best) {
			best <- s
			best_path <- p
		}
	}
	unname(best_path)
}

test_that("Viterbi.dthmm.adj returns the neutral state 3 for fewer than 2 genes", {
	expect_equal(Viterbi.dthmm.adj(hr_dthmm(c(0.5))), 3)
	expect_equal(Viterbi.dthmm.adj(hr_dthmm(numeric(0))), 3)
})

test_that("Viterbi.dthmm.adj follows a long run of values at one state's mean", {
	for (s in 1:6) {
		x <- rep(hr_means6[s], 12)
		expect_equal(Viterbi.dthmm.adj(hr_dthmm(x)), rep(s, 12))
	}
})

test_that("Viterbi.dthmm.adj finds a loss block between neutral flanks", {
	## 12 neutral, 12 at the 0.5 mean (state 2), 12 neutral; a switch costs
	## log(1e-6) per transition and each gene gains about 3 in log score
	x <- c(rep(1, 12), rep(0.5, 12), rep(1, 12))
	expect_equal(Viterbi.dthmm.adj(hr_dthmm(x)), rep(c(3, 2, 3), each=12))
	## gain block at 1.5 (state 4)
	x <- c(rep(1, 12), rep(1.5, 12), rep(1, 12))
	expect_equal(Viterbi.dthmm.adj(hr_dthmm(x)), rep(c(3, 4, 3), each=12))
})

test_that("Viterbi.dthmm.adj does not switch for a block too short to pay for two transitions", {
	## 2 genes of evidence cannot beat two transitions of cost log(1e-6) each
	x <- c(rep(1, 12), rep(0.5, 2), rep(1, 12))
	expect_equal(Viterbi.dthmm.adj(hr_dthmm(x)), rep(3, 26))
})

test_that("Viterbi.dthmm.adj matches an exhaustive search on short sequences", {
	## the first transition probability is large here so that paths with several
	## switches compete; checks the DP, the traceback and the emission scoring
	t <- 0.15
	h <- .get_HMM(hr_cnv_mean_sd(0.3), t)
	cases <- list(c(0.2, 0.9, 1.6, 1.0, 0.4),
				  c(1.0, 1.0, 0.4, 0.5, 2.9),
				  c(2.1, 0.1, 0.0, 1.4, 1.1),
				  c(0.7, 1.2))
	for (x in cases) {
		obj <- HiddenMarkov::dthmm(x=x, Pi=h$state_transitions, delta=h$delta, distn="norm",
								   pm=h$state_emission_params)
		sdm <- median(h$state_emission_params$sd)
		expect_equal(Viterbi.dthmm.adj(obj),
					 hr_bruteforce_viterbi(x, hr_means6, sdm, h$state_transitions, h$delta))
	}
})

test_that("Viterbi.dthmm.adj uses the median of the state sds for all states", {
	h <- .get_HMM(hr_cnv_mean_sd(0.1), 1e-6)
	x <- c(rep(1, 12), rep(0.5, 12), rep(1, 12))
	pm_a <- h$state_emission_params
	pm_b <- h$state_emission_params
	## same median (0.1) but one state with a huge sd: identical trace
	pm_b$sd[1] <- 50
	oa <- HiddenMarkov::dthmm(x=x, Pi=h$state_transitions, delta=h$delta, distn="norm", pm=pm_a)
	ob <- HiddenMarkov::dthmm(x=x, Pi=h$state_transitions, delta=h$delta, distn="norm", pm=pm_b)
	expect_equal(Viterbi.dthmm.adj(oa), Viterbi.dthmm.adj(ob))
})

## ---------------------------------------------------------------------------
## i6 prediction wrappers
## ---------------------------------------------------------------------------

## expected states for the group-A cells: loss block (state 2) on genes 13-24
hr_expected_loss_col <- function() {
	c(rep(3, 12), rep(2, 12), rep(3, 12), rep(3, 12))
}

test_that("predict_CNV_via_HMM_on_indiv_cells calls each cell and chromosome separately", {
	obj <- hr_make_obj()
	out <- predict_CNV_via_HMM_on_indiv_cells(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1))
	expect_equal(dim(out@expr.data), dim(obj@expr.data))
	expect_equal(dimnames(out@expr.data), dimnames(obj@expr.data))
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(3, 48))
	}
	## the input object is not modified
	expect_equal(obj@expr.data[13, 5], 0.5)
})

test_that("predict_CNV_via_HMM_on_indiv_cells assigns a state to every entry", {
	obj <- hr_make_obj()
	out <- predict_CNV_via_HMM_on_indiv_cells(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1))
	expect_true(all(out@expr.data %in% 1:6))
})

test_that("predict_CNV_via_HMM_on_indiv_cells treats a one-gene chromosome as neutral state 3", {
	obj <- hr_make_obj()
	## make c2 a single gene by moving 11 of its genes to a second one-gene-each layout
	obj@gene_order$chr[37:47] <- paste0("solo", 1:11)
	obj@expr.data[37:48, ] <- 0.5
	out <- predict_CNV_via_HMM_on_indiv_cells(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1))
	## every one-gene chromosome gets 3 whatever its value (Viterbi.dthmm.adj, n < 2)
	expect_equal(unname(out@expr.data[37:47, 1]), rep(3, 11))
})

test_that("predict_CNV_via_HMM_on_whole_tumor_samples runs per group when cluster_by_groups is TRUE", {
	obj <- hr_make_obj()
	out <- predict_CNV_via_HMM_on_whole_tumor_samples(obj, cluster_by_groups=TRUE,
													  cnv_mean_sd=hr_cnv_mean_sd(0.1),
													  cnv_level_to_mean_sd_fit=hr_sd_fits())
	## 4-cell groups: sd = 0.2/sqrt(4) = 0.1
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(3, 48))
	}
	expect_equal(dimnames(out@expr.data), dimnames(obj@expr.data))
})

test_that("predict_CNV_via_HMM_on_whole_tumor_samples averages the cells of a group before calling", {
	obj <- hr_make_obj()
	## one A cell sits at 0.6; the group mean over genes 13-24 is (3*0.5+0.6)/4 = 0.525,
	## 0.025 from the state-2 mean and 0.475 from neutral (sd 0.1 for 4 cells)
	obj@expr.data[13:24, 8] <- 0.6
	out <- predict_CNV_via_HMM_on_whole_tumor_samples(obj, TRUE, hr_cnv_mean_sd(0.1), hr_sd_fits())
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
})

test_that("predict_CNV_via_HMM_on_tumor_subclusters calls each subcluster separately", {
	sc <- list(subclusters=list(ref=list(ref_s1=1:4),
								A=list(A_s1=5:8),
								B=list(B_s1=9:12)))
	obj <- hr_make_obj(subclusters=sc)
	out <- predict_CNV_via_HMM_on_tumor_subclusters(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1),
													cnv_level_to_mean_sd_fit=hr_sd_fits())
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(3, 48))
	}
})

test_that("predict_CNV_via_HMM_on_tumor_subclusters splits a group into subclusters", {
	## A split in two subclusters of 2 cells (sd 0.2/sqrt(2)): only the
	## first has the loss, the second has been made neutral
	sc <- list(subclusters=list(ref=list(ref_s1=1:4),
								A=list(A_s1=5:6, A_s2=7:8),
								B=list(B_s1=9:12)))
	obj <- hr_make_obj(subclusters=sc)
	obj@expr.data[13:24, 7:8] <- 1
	out <- predict_CNV_via_HMM_on_tumor_subclusters(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1),
													cnv_level_to_mean_sd_fit=hr_sd_fits())
	for (cell in 5:6) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
	for (cell in 7:8) {
		expect_equal(unname(out@expr.data[, cell]), rep(3, 48))
	}
})

test_that("predict_CNV_via_HMM_on_tumor_subclusters_per_chr uses per-chromosome subclusters", {
	sc <- list(subclusters=list(ref=list(ref_s1=1:4), A=list(A_s1=5:8), B=list(B_s1=9:12)))
	obj <- hr_make_obj(subclusters=sc)
	## on c1 the cells are split as 5:6 / 7:8 / 9:12 / ref, on c2 as the whole objects
	per_chr <- list(c1=list(r=1:4, a1=5:6, a2=7:8, b=9:12),
					c2=list(r=1:4, a=5:8, b=9:12))
	out <- predict_CNV_via_HMM_on_tumor_subclusters_per_chr(obj, per_chr,
															cnv_mean_sd=hr_cnv_mean_sd(0.1),
															cnv_level_to_mean_sd_fit=hr_sd_fits())
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(3, 48))
	}
})

## ---------------------------------------------------------------------------
## state conversion
## ---------------------------------------------------------------------------

test_that("assign_HMM_states_to_proxy_expr_vals maps the six i6 states to CNV levels", {
	m <- matrix(c(1, 2, 3, 4, 5, 6, 3, 3), nrow=2, byrow=TRUE,
				dimnames=list(c("a", "b"), paste0("c", 1:4)))
	obj <- new("infercnv", expr.data=m, count.data=m,
			   gene_order=data.frame(chr=c("x", "x"), start=1:2, stop=1:2),
			   reference_grouped_cell_indices=list(), observation_grouped_cell_indices=list(),
			   tumor_subclusters=NULL, options=list(), .hspike=NULL)
	out <- assign_HMM_states_to_proxy_expr_vals(obj)
	expect_equal(out@expr.data, matrix(c(0, 0.5, 1, 1.5, 2, 3, 1, 1), nrow=2, byrow=TRUE,
									   dimnames=dimnames(m)))
})

test_that("assign_HMM_states_to_proxy_expr_vals does not remap an already converted value twice", {
	## 1 -> 0 and 3 -> 1: the new 1 (from state 3) must stay 1, not become 0
	m <- matrix(c(3, 1, 3, 3), nrow=1, dimnames=list("a", paste0("c", 1:4)))
	obj <- new("infercnv", expr.data=m, count.data=m,
			   gene_order=data.frame(chr="x", start=1, stop=1),
			   reference_grouped_cell_indices=list(), observation_grouped_cell_indices=list(),
			   tumor_subclusters=NULL, options=list(), .hspike=NULL)
	out <- assign_HMM_states_to_proxy_expr_vals(obj)
	expect_equal(unname(out@expr.data[1, ]), c(1, 0, 1, 1))
})

test_that("i3HMM_assign_HMM_states_to_proxy_expr_vals maps the three i3 states", {
	m <- matrix(c(1, 2, 3, 2, 1, 3), nrow=2, byrow=TRUE,
				dimnames=list(c("a", "b"), paste0("c", 1:3)))
	obj <- new("infercnv", expr.data=m, count.data=m,
			   gene_order=data.frame(chr=c("x", "x"), start=1:2, stop=1:2),
			   reference_grouped_cell_indices=list(), observation_grouped_cell_indices=list(),
			   tumor_subclusters=NULL, options=list(), .hspike=NULL)
	out <- i3HMM_assign_HMM_states_to_proxy_expr_vals(obj)
	expect_equal(out@expr.data, matrix(c(0.5, 1, 1.5, 1, 0.5, 1.5), nrow=2, byrow=TRUE,
									   dimnames=dimnames(m)))
})

## ---------------------------------------------------------------------------
## consensus and region definition
## ---------------------------------------------------------------------------

test_that(".get_state_consensus returns the most frequent state per gene as numeric", {
	m <- rbind(g1=c(3, 3, 2),     # 3
			   g2=c(2, 2, 2),     # 2
			   g3=c(1, 4, 4),     # 4
			   g4=c(5, 5, 6),     # 5
			   g5=c(6, 6, 6))     # 6
	out <- .get_state_consensus(m)
	expect_true(is.numeric(out))
	expect_equal(out, c(3, 2, 4, 5, 6))
})

test_that(".get_state_consensus works for a single cell and for a single gene", {
	expect_equal(.get_state_consensus(matrix(c(3, 2, 1), ncol=1)), c(3, 2, 1))
	expect_equal(.get_state_consensus(matrix(c(4, 4, 2), nrow=1)), 4)
})

test_that(".define_cnv_gene_regions splits a chromosome at each state change", {
	ge <- data.frame(chr=rep("chr1", 7), start=c(10, 20, 30, 40, 50, 60, 70),
					 stop=c(15, 25, 35, 45, 55, 65, 75), stringsAsFactors=FALSE)
	rownames(ge) <- paste0("g", 1:7)
	cons <- c(3, 3, 2, 2, 2, 3, 4)
	regs <- .define_cnv_gene_regions(cons, ge, 0)
	expect_equal(names(regs), paste0("chr1-region_", 1:4))
	expect_equal(regs[[1]]$gene, c("g1", "g2"))
	expect_equal(regs[[1]]$state, c(3, 3))
	expect_equal(regs[[1]]$start, c(10, 20))
	expect_equal(regs[[1]]$end, c(15, 25))
	expect_equal(regs[[1]]$chr, rep("chr1", 2))
	expect_equal(regs[[2]]$gene, c("g3", "g4", "g5"))
	expect_equal(regs[[2]]$state, c(2, 2, 2))
	expect_equal(regs[[3]]$gene, "g6")
	expect_equal(regs[[4]]$gene, "g7")
	expect_equal(regs[[4]]$state, 4)
	expect_equal(colnames(regs[[1]]), c("state", "gene", "chr", "start", "end"))
})

test_that(".define_cnv_gene_regions numbers regions continuously across chromosomes from the counter", {
	ge <- data.frame(chr=c("a", "a", "a", "b", "b"), start=1:5, stop=1:5 + 1, stringsAsFactors=FALSE)
	rownames(ge) <- paste0("g", 1:5)
	cons <- c(3, 2, 2, 3, 3)
	regs <- .define_cnv_gene_regions(cons, ge, 10)
	## chr a: region 11 (g1), region 12 (g2,g3); chr b: region 13 (g4,g5)
	expect_equal(names(regs), c("a-region_11", "a-region_12", "b-region_13"))
	expect_equal(regs[["b-region_13"]]$gene, c("g4", "g5"))
})

test_that(".define_cnv_gene_regions skips chromosomes with fewer than two genes", {
	ge <- data.frame(chr=c("a", "a", "solo", "b", "b"), start=1:5, stop=1:5, stringsAsFactors=FALSE)
	rownames(ge) <- paste0("g", 1:5)
	regs <- .define_cnv_gene_regions(c(3, 3, 5, 2, 2), ge, 0)
	## no region and no counter increment for the single-gene chromosome
	expect_equal(names(regs), c("a-region_1", "b-region_2"))
	expect_false(any(sapply(regs, function(r) "g3" %in% r$gene)))
})

test_that(".define_cnv_gene_regions keeps a state change at a chromosome boundary as separate regions", {
	ge <- data.frame(chr=c("a", "a", "b", "b"), start=c(1, 2, 1, 2), stop=c(1, 2, 1, 2),
					 stringsAsFactors=FALSE)
	rownames(ge) <- paste0("g", 1:4)
	## same state on both sides of the boundary: still one region per chromosome
	regs <- .define_cnv_gene_regions(c(3, 3, 3, 3), ge, 0)
	expect_equal(names(regs), c("a-region_1", "b-region_2"))
})

test_that(".define_cnv_gene_regions returns an empty list when no chromosome has two genes", {
	ge <- data.frame(chr=c("a", "b"), start=1:2, stop=1:2, stringsAsFactors=FALSE)
	rownames(ge) <- c("g1", "g2")
	expect_equal(.define_cnv_gene_regions(c(3, 2), ge, 0), list())
})

test_that(".get_cnv_gene_region_bounds reports min start and max end per region", {
	regs <- list("chr1-region_1"=data.frame(state=c(3, 3), gene=c("a", "b"), chr="chr1",
											start=c(100, 20), end=c(150, 60)),
				 "chr1-region_2"=data.frame(state=2, gene="c", chr="chr1", start=300, end=350))
	b <- .get_cnv_gene_region_bounds(regs)
	expect_equal(colnames(b), c("cnv_name", "state", "chr", "start", "end"))
	expect_equal(b$cnv_name, c("chr1-region_1", "chr1-region_2"))
	expect_equal(b$state, c(3, 2))
	expect_equal(b$start, c(20, 300))
	expect_equal(b$end, c(150, 350))
	expect_equal(rownames(b), c("1", "2"))
})

test_that(".get_cnv_gene_region_bounds on an empty region list returns NULL", {
	expect_null(.get_cnv_gene_region_bounds(list()))
})

## ---------------------------------------------------------------------------
## get_predicted_CNV_regions
## ---------------------------------------------------------------------------

test_that("get_predicted_CNV_regions by consensus calls each reference and observation group", {
	obj <- hr_make_state_obj()
	regs <- get_predicted_CNV_regions(obj, by="consensus")
	expect_length(regs, 2)
	expect_equal(sapply(regs, function(r) r$cell_group_name), c("ref", "tumA"))
	expect_equal(regs[[1]]$cells, c("ref1", "ref2"))
	expect_equal(regs[[2]]$cells, c("tumA1", "tumA2", "tumA3"))

	## ref group: all neutral; chr1, chr2 have regions, chr3 (one gene) is skipped
	expect_equal(names(regs[[1]]$gene_regions), c("chr1-region_1", "chr2-region_2"))
	expect_equal(regs[[1]]$cnv_ranges,
				 data.frame(cnv_name=c("chr1-region_1", "chr2-region_2"), state=c(3, 3),
							chr=c("chr1", "chr2"), start=c(100, 10), end=c(550, 35),
							stringsAsFactors=FALSE))

	## tumor group, per-gene mode over 3 cells:
	## chr1 genes 1-5: 3,2,2,3,3 (cell 3 supports 3 only where it differs: g5 has 3,3,4 -> 3)
	## chr2 genes 6-8: 5,5,3 (3,3,3 minority cell tumA3 loses to two cells)
	tum <- regs[[2]]
	expect_equal(names(tum$gene_regions),
				 c("chr1-region_3", "chr1-region_4", "chr1-region_5", "chr2-region_6", "chr2-region_7"))
	expect_equal(tum$cnv_ranges$state, c(3, 2, 3, 5, 3))
	expect_equal(tum$cnv_ranges$start, c(100, 200, 400, 10, 30))
	expect_equal(tum$cnv_ranges$end, c(150, 350, 550, 25, 35))
	expect_equal(tum$cnv_ranges$chr, c("chr1", "chr1", "chr1", "chr2", "chr2"))
})

test_that("get_predicted_CNV_regions by cell makes one entry per cell, named by column", {
	obj <- hr_make_state_obj()
	regs <- get_predicted_CNV_regions(obj, by="cell")
	## no subclusters defined: by falls back to consensus
	expect_equal(sapply(regs, function(r) r$cell_group_name), c("ref", "tumA"))

	sc <- list(subclusters=list(ref=list(ref_s1=1:2), tumA=list(tumA_s1=3:4, tumA_s2=5)))
	obj <- hr_make_state_obj(subclusters=sc)
	regs <- get_predicted_CNV_regions(obj, by="cell")
	expect_equal(sapply(regs, function(r) r$cell_group_name),
				 c("ref1", "ref2", "tumA1", "tumA2", "tumA3"))
	expect_equal(regs[[1]]$cells, "ref1")
	## cell tumA3 alone: chr1 3,3,3,3,4 -> regions {g1-4: 3}, {g5: 4}; chr2 all 3; counter
	## continues from the earlier cells: ref1, ref2 have 2 regions each, tumA1 and tumA2
	## have 5 each (chr1: 3,2,2,3,3 -> 3 regions; chr2: 5,5,3 -> 2 regions) = 14
	r5 <- regs[[5]]
	expect_equal(names(r5$gene_regions),
				 c("chr1-region_15", "chr1-region_16", "chr2-region_17"))
	expect_equal(r5$cnv_ranges$state, c(3, 4, 3))
	expect_equal(r5$cnv_ranges$start, c(100, 500, 10))
	expect_equal(r5$cnv_ranges$end, c(450, 550, 35))
})

test_that("get_predicted_CNV_regions by subcluster uses the tumor_subclusters definition", {
	sc <- list(subclusters=list(ref=list(ref_s1=1:2), tumA=list(tumA_s1=3:4, tumA_s2=5)))
	obj <- hr_make_state_obj(subclusters=sc)
	regs <- get_predicted_CNV_regions(obj, by="subcluster")
	expect_equal(sapply(regs, function(r) r$cell_group_name), c("ref.ref_s1", "tumA.tumA_s1", "tumA.tumA_s2"))
	expect_equal(regs[[2]]$cells, c("tumA1", "tumA2"))
	expect_equal(regs[[3]]$cells, "tumA3")
	## s1 (cells identical) keeps the loss block on genes 2-3 and the gain at chr2 genes 6-7
	expect_equal(regs[[2]]$cnv_ranges$state, c(3, 2, 3, 5, 3))
	expect_equal(regs[[2]]$cnv_ranges$start, c(100, 200, 400, 10, 30))
	## region names are unique across subclusters: counter continues
	all_names <- unlist(lapply(regs, function(r) names(r$gene_regions)))
	expect_equal(anyDuplicated(all_names), 0L)
	expect_equal(names(regs[[1]]$gene_regions), c("chr1-region_1", "chr2-region_2"))
	expect_equal(names(regs[[2]]$gene_regions)[1], "chr1-region_3")
})

test_that("get_predicted_CNV_regions rejects an unknown by value", {
	expect_error(get_predicted_CNV_regions(hr_make_state_obj(), by="bogus"))
})

## ---------------------------------------------------------------------------
## generate_cnv_region_reports
## ---------------------------------------------------------------------------

hr_read <- function(dir, prefix, ext) {
	read.table(file.path(dir, paste0(prefix, ext)), header=TRUE, sep="\t", stringsAsFactors=FALSE,
			   comment.char="")
}

test_that("generate_cnv_region_reports writes the four report files with the expected contents", {
	out_dir <- tempfile("hr_reports_")
	dir.create(out_dir)
	on.exit(unlink(out_dir, recursive=TRUE), add=TRUE)
	obj <- hr_make_state_obj()
	generate_cnv_region_reports(obj, "pre", out_dir, by="consensus")

	expect_setequal(list.files(out_dir),
					c("pre.cell_groupings", "pre.pred_cnv_regions.dat", "pre.pred_cnv_genes.dat",
					  "pre.genes_used.dat"))

	groups <- hr_read(out_dir, "pre", ".cell_groupings")
	expect_equal(colnames(groups), c("cell_group_name", "cell"))
	expect_equal(groups$cell_group_name, c("ref", "ref", rep("tumA", 3)))
	expect_equal(groups$cell, c("ref1", "ref2", "tumA1", "tumA2", "tumA3"))

	regions <- hr_read(out_dir, "pre", ".pred_cnv_regions.dat")
	expect_equal(colnames(regions), c("cell_group_name", "cnv_name", "state", "chr", "start", "end"))
	expect_equal(regions$cell_group_name, c("ref", "ref", rep("tumA", 5)))
	expect_equal(regions$cnv_name,
				 c("chr1-region_1", "chr2-region_2", "chr1-region_3", "chr1-region_4",
				   "chr1-region_5", "chr2-region_6", "chr2-region_7"))
	expect_equal(regions$state, c(3, 3, 3, 2, 3, 5, 3))
	expect_equal(regions$start, c(100, 10, 100, 200, 400, 10, 30))
	expect_equal(regions$end, c(550, 35, 150, 350, 550, 25, 35))

	genes <- hr_read(out_dir, "pre", ".pred_cnv_genes.dat")
	expect_equal(colnames(genes), c("cell_group_name", "gene_region_name", "state", "gene", "chr",
									"start", "end"))
	## ref: 5 + 3 genes; tumA: 8 genes in chr1/chr2 (chr3 has a single gene and is skipped)
	expect_equal(nrow(genes), 16)
	expect_equal(genes$gene[genes$cell_group_name == "tumA" & genes$gene_region_name == "chr1-region_4"],
				 c("gene2", "gene3"))
	expect_equal(genes$state[genes$cell_group_name == "tumA" & genes$gene_region_name == "chr2-region_6"],
				 c(5, 5))

	used <- read.table(file.path(out_dir, "pre.genes_used.dat"), header=TRUE, sep="\t",
					   row.names=1, stringsAsFactors=FALSE)
	expect_equal(rownames(used), paste0("gene", 1:9))
	expect_equal(used$chr, obj@gene_order$chr)
	expect_equal(used$start, obj@gene_order$start)
	expect_equal(used$stop, obj@gene_order$stop)
})

test_that("generate_cnv_region_reports drops the neutral state when ignore_neutral_state is given", {
	out_dir <- tempfile("hr_reports_")
	dir.create(out_dir)
	on.exit(unlink(out_dir, recursive=TRUE), add=TRUE)
	generate_cnv_region_reports(hr_make_state_obj(), "pre", out_dir, ignore_neutral_state=3, by="consensus")

	regions <- hr_read(out_dir, "pre", ".pred_cnv_regions.dat")
	expect_equal(regions$cnv_name, c("chr1-region_4", "chr2-region_6"))
	expect_equal(regions$state, c(2, 5))
	expect_equal(regions$cell_group_name, c("tumA", "tumA"))

	genes <- hr_read(out_dir, "pre", ".pred_cnv_genes.dat")
	expect_true(all(genes$state != 3))
	expect_equal(genes$gene, c("gene2", "gene3", "gene6", "gene7"))
	expect_equal(genes$gene_region_name, c("chr1-region_4", "chr1-region_4", "chr2-region_6", "chr2-region_6"))

	## the cell groupings and the genes used are not filtered
	groups <- hr_read(out_dir, "pre", ".cell_groupings")
	expect_equal(nrow(groups), 5)
	used <- read.table(file.path(out_dir, "pre.genes_used.dat"), header=TRUE, sep="\t", row.names=1)
	expect_equal(nrow(used), 9)
})

test_that("generate_cnv_region_reports by subcluster writes one row group per subcluster", {
	out_dir <- tempfile("hr_reports_")
	dir.create(out_dir)
	on.exit(unlink(out_dir, recursive=TRUE), add=TRUE)
	sc <- list(subclusters=list(ref=list(ref_s1=1:2), tumA=list(tumA_s1=3:4, tumA_s2=5)))
	generate_cnv_region_reports(hr_make_state_obj(subclusters=sc), "sub", out_dir, by="subcluster")
	groups <- hr_read(out_dir, "sub", ".cell_groupings")
	expect_equal(groups$cell_group_name,
				 c("ref.ref_s1", "ref.ref_s1", "tumA.tumA_s1", "tumA.tumA_s1", "tumA.tumA_s2"))
	regions <- hr_read(out_dir, "sub", ".pred_cnv_regions.dat")
	expect_equal(unique(regions$cell_group_name), c("ref.ref_s1", "tumA.tumA_s1", "tumA.tumA_s2"))
	## tumA_s2 is the single cell tumA3: chr1 3,3,3,3,4 and chr2 all 3
	s2 <- regions[regions$cell_group_name == "tumA.tumA_s2", ]
	expect_equal(s2$state, c(3, 4, 3))
	expect_equal(s2$start, c(100, 500, 10))
})

test_that("adjust_genes_regions_report keeps surviving regions and rewrites their states", {
	out_dir <- tempfile("hr_reports_")
	dir.create(out_dir)
	on.exit(unlink(out_dir, recursive=TRUE), add=TRUE)
	generate_cnv_region_reports(hr_make_state_obj(), "pre", out_dir, ignore_neutral_state=3, by="consensus")

	## the bayesian step kept chr1-region_4 (state changed to 1) and dropped chr2-region_6
	bn <- new("MCMC_inferCNV", expr.data=matrix(0), count.data=matrix(0),
			  gene_order=data.frame(), reference_grouped_cell_indices=list(),
			  observation_grouped_cell_indices=list(), tumor_subclusters=NULL,
			  cell_gene=list(list(cnv_regions="chr1-region_4", State=1)))
	adjust_genes_regions_report(bn, "pre", "post", out_dir)

	genes <- hr_read(out_dir, "post", ".pred_cnv_genes.dat")
	expect_equal(genes$gene_region_name, c("chr1-region_4", "chr1-region_4"))
	expect_equal(genes$gene, c("gene2", "gene3"))
	expect_equal(genes$state, c(1, 1))

	regions <- hr_read(out_dir, "post", ".pred_cnv_regions.dat")
	expect_equal(regions$cnv_name, "chr1-region_4")
	expect_equal(regions$state, 1)
	expect_equal(regions$start, 200)
	expect_equal(regions$end, 350)
})

test_that("adjust_genes_regions_report stops when the input report is missing", {
	out_dir <- tempfile("hr_reports_")
	dir.create(out_dir)
	on.exit(unlink(out_dir, recursive=TRUE), add=TRUE)
	bn <- new("MCMC_inferCNV", expr.data=matrix(0), count.data=matrix(0),
			  gene_order=data.frame(), reference_grouped_cell_indices=list(),
			  observation_grouped_cell_indices=list(), tumor_subclusters=NULL,
			  cell_gene=list())
	expect_error(adjust_genes_regions_report(bn, "absent", "post", out_dir), "Cannot find")
})

## ---------------------------------------------------------------------------
## end to end: HMM then regions
## ---------------------------------------------------------------------------

test_that("HMM states then region calling recover the planted loss block", {
	obj <- hr_make_obj()
	hmm <- predict_CNV_via_HMM_on_indiv_cells(obj, cnv_mean_sd=hr_cnv_mean_sd(0.1))
	regs <- get_predicted_CNV_regions(hmm, by="consensus")
	## groups: ref, A, B (reference first, then observations in list order)
	expect_equal(sapply(regs, function(r) r$cell_group_name), c("ref", "A", "B"))
	a <- regs[[2]]$cnv_ranges
	## c1: neutral g1-12 (start 1, end 13), loss g13-24 (start 13, end 25), neutral g25-36;
	## c2: neutral. ends are start+1 in this fixture.
	expect_equal(a$state, c(3, 2, 3, 3))
	expect_equal(a$chr, c("c1", "c1", "c1", "c2"))
	expect_equal(a$start, c(1, 13, 25, 37))
	expect_equal(a$end, c(13, 25, 37, 49))
	## reference and group B are neutral over the whole of each chromosome
	expect_equal(regs[[1]]$cnv_ranges$state, c(3, 3))
	expect_equal(regs[[3]]$cnv_ranges$state, c(3, 3))
})

## ---------------------------------------------------------------------------
## i3 HMM
## ---------------------------------------------------------------------------

test_that("determine_mean_delta_via_Z is |qnorm(p)| * sigma", {
	## qnorm(0.05) = -1.6448536269514722, qnorm(0.025) = -1.959963984540054
	expect_equal(determine_mean_delta_via_Z(1, 0.05), 1.6448536269514722, tolerance=1e-12)
	expect_equal(determine_mean_delta_via_Z(0.3, 0.05), 0.3 * 1.6448536269514722, tolerance=1e-12)
	expect_equal(determine_mean_delta_via_Z(2, 0.025), 2 * 1.959963984540054, tolerance=1e-12)
	expect_equal(determine_mean_delta_via_Z(2, 0.5), 0)
	## the absolute value makes p > 0.5 positive too
	expect_equal(determine_mean_delta_via_Z(1, 0.975), 1.959963984540054, tolerance=1e-12)
})

test_that("get_HoneyBADGER_setGexpDev is reproducible and decreases with the number of cells", {
	set.seed(7)
	d40 <- get_HoneyBADGER_setGexpDev(gexp.sd=0.5, alpha=0.05, k_cells=40, n_iter=30)
	set.seed(7)
	d40b <- get_HoneyBADGER_setGexpDev(gexp.sd=0.5, alpha=0.05, k_cells=40, n_iter=30)
	expect_equal(d40, d40b)
	## 40 cells per sample resolve a shift of well under one sd at p = 0.05
	expect_true(d40 > 0 && d40 < 0.5)
	set.seed(7)
	d10 <- get_HoneyBADGER_setGexpDev(gexp.sd=0.5, alpha=0.05, k_cells=10, n_iter=30)
	expect_true(d10 > d40)
})

test_that("get_HoneyBADGER_setGexpDev scales linearly with the sd for a fixed seed", {
	## devs = seq(0, sd, sd/10) and the samples are rnorm(., 0, sd) / rnorm(., dev, sd): the
	## whole simulation is the sd=1 simulation multiplied by sd (same draws, KS is scale free)
	set.seed(11)
	a <- get_HoneyBADGER_setGexpDev(gexp.sd=1, alpha=0.05, k_cells=20, n_iter=20)
	set.seed(11)
	b <- get_HoneyBADGER_setGexpDev(gexp.sd=3, alpha=0.05, k_cells=20, n_iter=20)
	expect_equal(unname(b), unname(3 * a), tolerance=1e-8)
})

test_that("get_HoneyBADGER_setGexpDev treats k_cells below 2 as 2", {
	set.seed(3)
	a <- get_HoneyBADGER_setGexpDev(gexp.sd=0.5, alpha=0.05, k_cells=1, n_iter=10)
	set.seed(3)
	b <- get_HoneyBADGER_setGexpDev(gexp.sd=0.5, alpha=0.05, k_cells=2, n_iter=10)
	expect_equal(a, b)
})

test_that(".i3HMM_get_sd_trend_by_num_cells_fit uses the reference cells' mean and sd", {
	obj <- hr_make_obj()
	## reference cells (1-4) get distinct values, so mean and sd are checkable by hand
	obj@expr.data[, 1:4] <- matrix(rep(c(0.5, 1, 1.5, 1), each=48), nrow=48)
	vals <- c(0.5, 1, 1.5, 1)
	set.seed(1)
	st <- suppressMessages(.i3HMM_get_sd_trend_by_num_cells_fit(obj, i3_p_val=0.05))
	expect_equal(names(st), c("mu", "sigma", "mean_delta", "KS_delta"))
	expect_equal(st$mu, 1)
	## 48 genes x each value; sum sq dev per gene row = 0.25+0+0.25+0 = 0.5, N = 192
	expect_equal(st$sigma, sqrt(48 * 0.5 / (192 - 1)), tolerance=1e-12)
	expect_equal(st$mean_delta, 1.6448536269514722 * st$sigma, tolerance=1e-12)
	expect_length(st$KS_delta, 1)
	expect_true(is.finite(st$KS_delta))
})

test_that(".i3HMM_get_sd_trend_by_num_cells_fit falls back to observations when there is no reference", {
	obj <- hr_make_obj()
	obj@reference_grouped_cell_indices <- list()
	obj@observation_grouped_cell_indices <- list(A=5:8)
	## A has 12 of 48 genes at 0.5 in all 4 cells: mean = (36*1 + 12*0.5)/48 = 0.875
	set.seed(1)
	st <- suppressMessages(.i3HMM_get_sd_trend_by_num_cells_fit(obj, i3_p_val=0.1))
	expect_equal(st$mu, 0.875, tolerance=1e-12)
	## 12 genes at deviation 0.375 and 36 at 0.125, per cell column: 12*0.140625 + 36*0.015625
	expect_equal(st$sigma, sqrt((12 * 0.375^2 + 36 * 0.125^2) * 4 / (192 - 1)), tolerance=1e-12)
	expect_equal(st$mean_delta, 1.2815515655446004 * st$sigma, tolerance=1e-12)
})

test_that(".i3HMM_get_HMM builds the 3-state model with mu -/+ delta emissions", {
	t <- 0.02
	st <- list(mu=1.2, sigma=0.3, mean_delta=0.4, KS_delta=0.7)
	h <- .i3HMM_get_HMM(st, t, use_KS=FALSE)
	expected_pi <- matrix(t, 3, 3)
	diag(expected_pi) <- 1 - 5 * t
	expect_equal(h$state_transitions, expected_pi)
	expect_equal(h$delta, c(t, 1 - 5 * t, t))
	expect_equal(h$state_emission_params$mean, c(0.8, 1.2, 1.6))
	expect_equal(h$state_emission_params$sd, rep(0.3, 3))
	## the KS delta is used when requested
	hk <- .i3HMM_get_HMM(st, t, use_KS=TRUE)
	expect_equal(hk$state_emission_params$mean, c(0.5, 1.2, 1.9))
})

hr_sd_trend <- function() list(mu=1, sigma=0.1, mean_delta=0.5, KS_delta=0.5)

## expected i3 states (1 loss, 2 neutral, 3 gain) for the planted loss in group A
hr_expected_loss_col_i3 <- function() c(rep(2, 12), rep(1, 12), rep(2, 12), rep(2, 12))

test_that("i3HMM_predict_CNV_via_HMM_on_indiv_cells calls a loss block as state 1", {
	obj <- hr_make_obj()
	obj@expr.data[13:24, 9:12] <- 1.5  # a gain in group B
	out <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_indiv_cells(obj, sd_trend=hr_sd_trend()))
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col_i3())
	}
	for (cell in 9:12) {
		expect_equal(unname(out@expr.data[, cell]), c(rep(2, 12), rep(3, 12), rep(2, 24)))
	}
	for (cell in 1:4) {
		expect_equal(unname(out@expr.data[, cell]), rep(2, 48))
	}
	expect_equal(dimnames(out@expr.data), dimnames(obj@expr.data))
})

test_that("i3HMM_predict_CNV_via_HMM_on_indiv_cells computes sd_trend from the object by default", {
	obj <- hr_make_obj()
	## reference cells need spread, otherwise sigma is 0
	obj@expr.data[, 1:4] <- 1 + 0.05 * outer(1:48, 1:4, function(i, j) c(-1, 0, 1, 0)[(i + j) %% 4 + 1])
	set.seed(5)
	a <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_indiv_cells(obj, i3_p_val=0.05))
	set.seed(5)
	st <- suppressMessages(.i3HMM_get_sd_trend_by_num_cells_fit(obj, 0.05))
	b <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_indiv_cells(obj, i3_p_val=0.05, sd_trend=st))
	expect_equal(a@expr.data, b@expr.data)
	expect_true(all(a@expr.data %in% 1:3))
})

test_that("i3HMM_predict_CNV_via_HMM_on_indiv_cells honours use_KS", {
	obj <- hr_make_obj()
	## a shallow loss (0.9) is only called when the loss mean sits close to the neutral one:
	## KS_delta 0.1 puts the del mean at 0.9, mean_delta 0.5 puts it at 0.5
	obj@expr.data[13:24, 5:8] <- 0.9
	st <- list(mu=1, sigma=0.01, mean_delta=0.5, KS_delta=0.1)
	ks <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_indiv_cells(obj, sd_trend=st, use_KS=TRUE))
	z <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_indiv_cells(obj, sd_trend=st, use_KS=FALSE))
	expect_equal(unname(ks@expr.data[, 5]), c(rep(2, 12), rep(1, 12), rep(2, 24)))
	expect_equal(unname(z@expr.data[, 5]), rep(2, 48))
})

test_that("i3HMM_predict_CNV_via_HMM_on_tumor_subclusters calls each subcluster", {
	sc <- list(subclusters=list(ref=list(ref_s1=1:4), A=list(A_s1=5:8), B=list(B_s1=9:12)))
	obj <- hr_make_obj(subclusters=sc)
	out <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(obj, sd_trend=hr_sd_trend()))
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col_i3())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(2, 48))
	}
})

test_that("i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples runs per group", {
	obj <- hr_make_obj()
	out <- suppressMessages(i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples(
		obj, cluster_by_groups=TRUE, sd_trend=hr_sd_trend()))
	for (cell in 5:8) {
		expect_equal(unname(out@expr.data[, cell]), hr_expected_loss_col_i3())
	}
	for (cell in c(1:4, 9:12)) {
		expect_equal(unname(out@expr.data[, cell]), rep(2, 48))
	}
})
