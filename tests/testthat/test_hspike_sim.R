## Tests for the hidden spike (R/inferCNV_hidden_spike.R), the mean-variance
## simulation (R/inferCNV_meanVarSim.R), the pieces of R/inferCNV_simple_sim.R
## that the h-spike path calls, and R/SplatterScrape.R.
##
## Most simulator tests use inputs whose answer does not depend on the random
## stream: a spline fitted to a constant variance of 0 makes every draw equal
## round(mean), and a constant dropout spline of 0 or 1 keeps or removes every
## value. The remaining tests check properties that hold for any draw, with
## tolerances sized from the number of draws.

local({
    old <- futile.logger::flog.threshold()
    futile.logger::flog.threshold(futile.logger::ERROR)
    withr::defer(futile.logger::flog.threshold(old), envir=testthat::teardown_env())
})

## ---- helpers (prefixed hs_) ----

hs_make_obj <- function(mat, ref, obs) {
    ## ref and obs are named lists of column indices; ref may be empty
    gene_order <- data.frame(chr=rep("chr1", nrow(mat)),
                             start=seq_len(nrow(mat)),
                             stop=seq_len(nrow(mat)),
                             row.names=rownames(mat),
                             stringsAsFactors=FALSE)
    new(Class="infercnv",
        expr.data=mat,
        count.data=mat,
        gene_order=gene_order,
        reference_grouped_cell_indices=ref,
        observation_grouped_cell_indices=obs)
}

## constant-per-gene matrix: every cell of a gene has the same value, so the
## per-group variance is 0 for every gene and the fitted spline is 0 everywhere.
hs_const_matrix <- function(ref_vals, obs_vals, n_ref=5, n_obs=5) {
    ng <- length(ref_vals)
    mat <- cbind(matrix(rep(ref_vals, n_ref), nrow=ng),
                 matrix(rep(obs_vals, n_obs), nrow=ng))
    rownames(mat) <- paste0("g", seq_len(ng))
    colnames(mat) <- paste0("c", seq_len(ncol(mat)))
    mat
}

## spline through a constant: predict() returns that constant everywhere
hs_const_spline <- function(value, from=0.1, to=60) {
    x <- log(seq(from, to, length.out=10))
    smooth.spline(x, rep(value, length(x)))
}

hs_dropout_params <- function(p) {
    list(spline=hs_const_spline(p))
}

hs_quiet <- function(expr) {
    ## the simple/meanvar code prints with cat() when the logistic fit fails
    out <- NULL
    utils::capture.output(out <- force(expr))
    out
}

hs_cnv_by_chr <- c(chrA=1, chr_0=0.01, chr_B=1, chr_0pt5=0.5, chr_C=1, chr_1pt5=1.5,
                   chr_D=1, chr_2pt0=2, chr_E=1, chr_3pt0=3, chr_F=1)

## ---- .get_hspike_chr_info ----

test_that(".get_hspike_chr_info lays out eleven fake chromosomes with fixed copy numbers", {
    ci <- .get_hspike_chr_info(400, 10000)
    expect_length(ci, 11)
    expect_equal(vapply(ci, function(x) x$name, ""), names(hs_cnv_by_chr))
    expect_equal(vapply(ci, function(x) x$cnv, 0), unname(hs_cnv_by_chr))
    ## 10000 - 4000 = 6000 genes remain and are all placed on the last chromosome
    expect_equal(vapply(ci, function(x) x$ngenes, 0), c(rep(400, 10), 6000))
})

test_that(".get_hspike_chr_info gives the last chromosome at least num_genes_each genes", {
    ## 5000 - 4000 = 1000 remaining
    expect_equal(.get_hspike_chr_info(400, 5000)[[11]]$ngenes, 1000)
    ## exactly 4400 total leaves 400
    expect_equal(.get_hspike_chr_info(400, 4400)[[11]]$ngenes, 400)
    ## fewer genes than 10 * num_genes_each: floor of num_genes_each
    expect_equal(.get_hspike_chr_info(400, 30)[[11]]$ngenes, 400)
    expect_equal(.get_hspike_chr_info(400, 4000)[[11]]$ngenes, 400)
    ## the dummy call used by .get_gene_expr_by_cnv
    ci <- .get_hspike_chr_info(1, 1)
    expect_equal(vapply(ci, function(x) x$ngenes, 0), rep(1, 11))
    expect_equal(vapply(ci, function(x) x$name, ""), names(hs_cnv_by_chr))
})

## ---- logistic and mean-vs-p0 tables ----

test_that(".logistic_midpt_slope is the logistic function", {
    expect_equal(.logistic_midpt_slope(2, midpt=2, slope=3), 0.5)
    expect_equal(.logistic_midpt_slope(2 + log(3), midpt=2, slope=1), 0.75)
    expect_equal(.logistic_midpt_slope(2 + log(3), midpt=2, slope=-1), 0.25)
    expect_equal(.logistic_midpt_slope(c(-1, 0, 1), midpt=0, slope=2),
                 1 / (1 + exp(c(2, 0, -2))))
    expect_equal(.logistic_midpt_slope(1000, 0, 1), 1)
    expect_equal(.logistic_midpt_slope(-1000, 0, 1), 0)
})

test_that(".get_mean_vs_p0_from_matrix returns row means and the fraction of zero cells", {
    m <- matrix(c(0, 0, 0, 0,
                  1, 0, 3, 0,
                  2, 2, 2, 2,
                  0, 4, 0, 0),
                nrow=4, byrow=TRUE, dimnames=list(paste0("g", 1:4), paste0("c", 1:4)))
    tb <- .get_mean_vs_p0_from_matrix(m)
    expect_s3_class(tb, "data.frame")
    expect_equal(colnames(tb), c("m", "p0"))
    expect_equal(rownames(tb), paste0("g", 1:4))
    expect_equal(tb$m, c(0, 1, 2, 1))
    expect_equal(tb$p0, c(1, 0.5, 0, 0.75))
})

test_that(".get_mean_vs_p0_table_from_matrix defaults to one group of all cells", {
    m <- matrix(c(0, 2, 4, 0,
                  1, 1, 1, 1), nrow=2, byrow=TRUE, dimnames=list(c("a", "b"), NULL))
    tb <- .get_mean_vs_p0_table_from_matrix(m)
    expect_equal(nrow(tb), 2)
    expect_equal(tb$m, c(1.5, 1))
    expect_equal(tb$p0, c(0.5, 0))
    expect_equal(tb$group_name, c("allcells", "allcells"))
})

test_that(".get_mean_vs_p0_table_from_matrix stacks groups in list order", {
    m <- matrix(c(0, 0, 4, 8,
                  1, 3, 0, 0), nrow=2, byrow=TRUE, dimnames=list(c("a", "b"), NULL))
    tb <- .get_mean_vs_p0_table_from_matrix(m, list(second=3:4, first=1:2))
    expect_equal(nrow(tb), 4)
    expect_equal(tb$group_name, c("second", "second", "first", "first"))
    expect_equal(tb$m, c(6, 0, 0, 2))
    expect_equal(tb$p0, c(0, 1, 1, 0))
})

test_that(".get_mean_vs_p0_table orders observation groups before reference groups", {
    m <- matrix(c(0, 0, 2, 2, 4, 4,
                  1, 1, 0, 0, 6, 6), nrow=2, byrow=TRUE,
                dimnames=list(c("a", "b"), paste0("c", 1:6)))
    obj <- hs_make_obj(m, ref=list(refA=1:2), obs=list(obsB=3:4, obsC=5:6))
    tb <- .get_mean_vs_p0_table(obj)
    expect_equal(tb$group_name, c("obsB", "obsB", "obsC", "obsC", "refA", "refA"))
    expect_equal(tb$m, c(2, 0, 4, 6, 0, 1))
    expect_equal(tb$p0, c(0, 1, 0, 0, 1, 0))
})

test_that(".get_logistic_params recovers the parameters of a logistic curve in log(mean)", {
    x <- seq(-2, 4, length.out=40)
    p0 <- 1 / (1 + exp(1.5 * (x - 0.5))) + 0.01 * sin(7 * x)
    tb <- data.frame(m=exp(x), p0=p0)
    lp <- hs_quiet(.get_logistic_params(tb))
    expect_setequal(names(lp), c("midpt", "slope", "spline"))
    expect_equal(lp$midpt, 0.5, tolerance=0.1)
    expect_equal(lp$slope, -1.5, tolerance=0.1)
    expect_s3_class(lp$spline, "smooth.spline")
    expect_equal(predict(lp$spline, x)$y, p0, tolerance=0.05)
})

test_that(".get_logistic_params ignores genes with zero mean", {
    x <- seq(-2, 4, length.out=40)
    p0 <- 1 / (1 + exp(1.5 * (x - 0.5))) + 0.01 * sin(7 * x)
    tb <- data.frame(m=c(exp(x), 0, 0), p0=c(p0, 1, 1))
    lp <- hs_quiet(.get_logistic_params(tb))
    ## the two zero-mean rows would otherwise put -Inf into log(m)
    expect_equal(lp$midpt, 0.5, tolerance=0.1)
    expect_equal(lp$slope, -1.5, tolerance=0.1)
    expect_equal(predict(lp$spline, x)$y, p0, tolerance=0.05)
})

test_that(".get_logistic_params falls back to the spline alone when the logistic cannot be fit", {
    ## p0 is exactly 0 everywhere, which a logistic cannot reach
    tb <- data.frame(m=c(1, 2, 4, 8, 16, 32, 64, 128), p0=0)
    lp <- hs_quiet(.get_logistic_params(tb))
    expect_equal(names(lp), "spline")
    expect_equal(predict(lp$spline, log(c(2, 10, 50)))$y, c(0, 0, 0), tolerance=1e-6)
})

test_that(".mean_vs_p0_to_stats fits a logistic and a spline to log(mean + 1)", {
    x <- log(seq(0.2, 20, length.out=40))
    p0 <- 1 / (1 + exp(1.5 * (x - 0.5))) + 0.01 * sin(7 * x)
    tb <- data.frame(m=exp(x) - 1, p0=pmin(pmax(p0, 0), 1))
    r <- .mean_vs_p0_to_stats(tb)
    expect_length(r, 6)
    expect_equal(r$logistic_x, log(tb$m + 1), tolerance=1e-8)
    expect_equal(r$logistic_y, tb$p0, tolerance=0.05)
    expect_s3_class(r[[5]], "smooth.spline")
    expect_s3_class(r[[6]], "nls")
    expect_equal(unname(coef(r[[6]])["x0"]), 0.5, tolerance=0.2)
})

test_that(".estimate_common_dispersion is larger for overdispersed counts than for Poisson counts", {
    skip_if_not_installed("edgeR")
    set.seed(21)
    mu <- rep(c(5, 20, 60, 120), 50)
    pois <- matrix(rpois(200 * 40, mu), nrow=200)
    nb <- matrix(rnbinom(200 * 40, mu=mu, size=2), nrow=200)  # dispersion 0.5
    d_pois <- .estimate_common_dispersion(pois)
    d_nb <- .estimate_common_dispersion(nb)
    expect_true(is.numeric(d_pois) && length(d_pois) == 1 && d_pois >= 0)
    expect_lt(d_pois, 0.1)
    expect_equal(d_nb, 0.5, tolerance=0.2)
})

## ---- mean / variance tables ----

test_that(".get_mean_var_given_matrix computes per-gene mean and sample variance per group", {
    m <- matrix(c(1, 3, 5, 7,
                  2, 2, 2, 2,
                  0, 4, 0, 8), nrow=3, byrow=TRUE,
                dimnames=list(c("a", "b", "c"), paste0("c", 1:4)))
    tb <- .get_mean_var_given_matrix(m)
    expect_s3_class(tb, "data.frame")
    expect_equal(colnames(tb), c("g", "m", "v"))
    expect_equal(tb$g, rep("allcells", 3))
    expect_equal(tb$m, c(4, 2, 3))
    ## var of 1,3,5,7 is 20/3; of 0,4,0,8 is 44/3
    expect_equal(tb$v, c(20 / 3, 0, 44 / 3))
    expect_equal(rownames(tb), c("a", "b", "c"))
})

test_that(".get_mean_var_given_matrix stacks groups in list order", {
    m <- matrix(c(1, 3, 10, 30,
                  2, 2, 5, 9), nrow=2, byrow=TRUE,
                dimnames=list(c("a", "b"), NULL))
    tb <- .get_mean_var_given_matrix(m, list(hi=3:4, lo=1:2))
    expect_equal(tb$g, c("hi", "hi", "lo", "lo"))
    expect_equal(tb$m, c(20, 7, 2, 2))
    expect_equal(tb$v, c(200, 8, 2, 0))
})

test_that(".get_mean_var_table orders observation groups before reference groups", {
    m <- matrix(c(1, 3, 10, 30, 0, 6,
                  2, 2, 5, 9, 4, 4), nrow=2, byrow=TRUE,
                dimnames=list(c("a", "b"), paste0("c", 1:6)))
    obj <- hs_make_obj(m, ref=list(refR=5:6), obs=list(obsX=1:2, obsY=3:4))
    tb <- .get_mean_var_table(obj)
    expect_equal(tb$g, c("obsX", "obsX", "obsY", "obsY", "refR", "refR"))
    expect_equal(tb$m, c(2, 2, 20, 7, 3, 4))
    expect_equal(tb$v, c(2, 0, 200, 8, 18, 0))
})

## ---- per-value simulation ----

test_that(".sim_expr_val_mean_var_no_dropout returns round(mean) when the fitted variance is 0", {
    sp <- hs_const_spline(0)
    expect_equal(.sim_expr_val_mean_var_no_dropout(3.2, sp), 3)
    expect_equal(.sim_expr_val_mean_var_no_dropout(10.7, sp), 11)
    expect_equal(.sim_expr_val_mean_var_no_dropout(0.2, sp), 0)
    expect_equal(.sim_expr_val_mean_var_no_dropout(250, sp), 250)
})

test_that(".sim_expr_val_mean_var_no_dropout treats a negative fitted log-variance as variance 0", {
    ## exp(-0.5) - 1 < 0 is clamped to 0
    sp <- hs_const_spline(-0.5)
    expect_equal(.sim_expr_val_mean_var_no_dropout(6.4, sp), 6)
})

test_that(".sim_expr_val_mean_var_no_dropout returns 0 without drawing for non-positive means", {
    sp <- hs_const_spline(0)
    set.seed(1)
    before <- .Random.seed
    expect_identical(.sim_expr_val_mean_var_no_dropout(0, sp), 0)
    expect_identical(.sim_expr_val_mean_var_no_dropout(-4, sp), 0)
    expect_identical(.Random.seed, before)
})

test_that(".sim_expr_val_mean_var_no_dropout draws round(max(rnorm(m, sqrt(var)), 0))", {
    ## spline of log(v + 1) = log(m + 1) is the identity, so the modelled variance equals the mean
    lm <- log(c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100) + 1)
    sp <- smooth.spline(lm, lm)
    for (m in c(0.3, 4, 37.5)) {
        v <- max(exp(predict(sp, log(m + 1))$y) - 1, 0)
        expect_equal(v, m, tolerance=1e-6)
        set.seed(99)
        expected <- round(max(rnorm(1, mean=m, sd=sqrt(v)), 0))
        set.seed(99)
        expect_identical(.sim_expr_val_mean_var_no_dropout(m, sp), expected)
    }
})

test_that(".sim_expr_val_mean_var_no_dropout draws non-negative integers with the modelled mean and variance", {
    lm <- log(c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100) + 1)
    sp <- smooth.spline(lm, lm)
    set.seed(5)
    draws <- vapply(seq_len(3000), function(i) .sim_expr_val_mean_var_no_dropout(50, sp), 0)
    expect_true(all(draws >= 0 & draws == round(draws)))
    ## sd of the mean of 3000 draws is sqrt(50 / 3000) = 0.13
    expect_equal(mean(draws), 50, tolerance=0.02)
    ## relative error of a variance estimate from 3000 draws is about 2.6%
    expect_equal(var(draws), 50, tolerance=0.15)
    ## a mean near 0 gives many zeros, never negatives
    small <- vapply(seq_len(500), function(i) .sim_expr_val_mean_var_no_dropout(0.05, sp), 0)
    expect_true(all(small >= 0))
    expect_gt(mean(small == 0), 0.9)
})

test_that(".sim_expr_val_mean_var without dropout parameters matches the no-dropout function", {
    lm <- log(c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100) + 1)
    sp <- smooth.spline(lm, lm)
    set.seed(8)
    a <- vapply(1:50, function(i) .sim_expr_val_mean_var(12, sp, NULL), 0)
    set.seed(8)
    b <- vapply(1:50, function(i) .sim_expr_val_mean_var_no_dropout(12, sp), 0)
    expect_identical(a, b)
    expect_identical(.sim_expr_val_mean_var(0, sp, NULL), 0)
})

test_that(".sim_expr_val_mean_var applies a dropout probability that is read from the spline at log(value)", {
    sp0 <- hs_const_spline(0)
    ## probability 1 removes every positive value
    expect_equal(.sim_expr_val_mean_var(5, sp0, hs_dropout_params(1)), 0)
    ## probability 0 keeps it
    expect_equal(.sim_expr_val_mean_var(5, sp0, hs_dropout_params(0)), 5)
    ## a mean of 0 stays 0 whatever the dropout parameters
    expect_equal(.sim_expr_val_mean_var(0, sp0, hs_dropout_params(0)), 0)
    set.seed(3)
    d <- vapply(1:4000, function(i) .sim_expr_val_mean_var(5, sp0, hs_dropout_params(0.5)), 0)
    expect_true(all(d %in% c(0, 5)))
    ## standard deviation of the zero fraction is 0.008
    expect_equal(mean(d == 0), 0.5, tolerance=0.1)
})

## ---- .apply_dropout ----

hs_count_matrix <- function(ng=30, nc=10, lambda=3, seed=1) {
    set.seed(seed)
    x <- matrix(rpois(ng * nc, lambda), ng, nc,
                dimnames=list(paste0("g", seq_len(ng)), paste0("c", seq_len(nc))))
    x[x == 0] <- 1  # every row has a positive mean and no zeros
    x
}

test_that(".apply_dropout with probability 0 returns the matrix unchanged, with dimnames", {
    x <- hs_count_matrix()
    r <- .apply_dropout(x, hs_dropout_params(0))
    expect_identical(dim(r), dim(x))
    expect_equal(r, x)
    expect_identical(dimnames(r), dimnames(x))
})

test_that(".apply_dropout with probability 1 zeroes every cell of every row", {
    x <- hs_count_matrix()
    r <- .apply_dropout(x, hs_dropout_params(1))
    expect_identical(dim(r), dim(x))
    expect_true(all(r == 0))
    expect_identical(dimnames(r), dimnames(x))
})

test_that(".apply_dropout leaves a row alone when it already has more zeros than the target", {
    ## target dropout 0.2 of 10 cells = 2 zeros; the row already has 5
    x <- hs_count_matrix(ng=8, nc=10)
    x[3, 1:5] <- 0
    set.seed(4)
    r <- .apply_dropout(x, hs_dropout_params(0.2))
    expect_identical(r[3, ], x[3, ])
})

test_that(".apply_dropout only turns values into zeros and reaches the target zero fraction", {
    ## dropout probability falls linearly with log(mean): p = 0.8 - 0.2 * log(mean)
    set.seed(12)
    lam <- seq(1, 6, length.out=300)
    nc <- 100
    x <- t(vapply(lam, function(l) rpois(nc, l), numeric(nc)))
    rownames(x) <- paste0("g", seq_len(nrow(x)))
    keep_rows <- rowSums(x) > 0
    x <- x[keep_rows, ]
    lm <- seq(-1, 3, length.out=30)
    dsp <- list(spline=smooth.spline(lm, 0.8 - 0.2 * lm))
    set.seed(13)
    r <- .apply_dropout(x, dsp)
    expect_identical(dim(r), dim(x))
    expect_true(all(r == 0 | r == x))
    expect_true(all(r[x == 0] == 0))
    ## reference expectation: row i ends with max(p_i * nc, zeros_i) zeros on average
    p <- 0.8 - 0.2 * log(rowMeans(x))
    nz <- rowSums(x == 0)
    expected_frac <- mean(pmax(p * nc, nz) / nc)
    ## sd of the overall zero fraction is below 0.005
    expect_equal(mean(r == 0), expected_frac, tolerance=0.03)
    ## rows that had no zeros and whose probability is high lose a large share of their values
    expect_gt(mean(r[which(nz == 0 & p > 0.5), ] == 0), 0.4)
})

## ---- meanvar matrix simulation ----

test_that(".get_simulated_cell_matrix_using_meanvar_trend_helper simulates round(mean) when the trend variance is 0", {
    m <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50)
    mvt <- data.frame(g="x", m=m, v=0)
    gm <- c(a=3.2, b=0, c=10.7, d=0.2)
    r <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gm, mvt, 3)
    expect_true(is.matrix(r))
    expect_identical(dim(r), c(4L, 3L))
    expect_identical(rownames(r), c("a", "b", "c", "d"))
    expect_identical(colnames(r), c("sim_cell_1", "sim_cell_2", "sim_cell_3"))
    expect_equal(r[, 1], c(a=3, b=0, c=11, d=0))
    expect_equal(r[, 2], r[, 1])
    expect_equal(r[, 3], r[, 1])
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend_helper applies dropout parameters when given", {
    m <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50)
    mvt <- data.frame(g="x", m=m, v=0)
    gm <- c(a=3.2, b=7.9, c=10.7, d=2.2)
    all_gone <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gm, mvt, 6, hs_dropout_params(1))
    expect_identical(dim(all_gone), c(4L, 6L))
    expect_true(all(all_gone == 0))
    none_lost <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gm, mvt, 6, hs_dropout_params(0))
    expect_equal(none_lost[, 1], c(a=3, b=8, c=11, d=2))
    expect_equal(none_lost, .get_simulated_cell_matrix_using_meanvar_trend_helper(gm, mvt, 6))
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend_helper follows the fitted mean-variance trend", {
    ## v = m, so log(v + 1) = log(m + 1) is fit by the identity spline
    m <- c(0.1, 0.5, 1, 2, 5, 10, 20, 50, 100)
    mvt <- data.frame(g="x", m=m, v=m)
    gm <- c(lo=4, hi=50)
    set.seed(31)
    r <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gm, mvt, 2000)
    expect_true(all(r >= 0 & r == round(r)))
    ## standard error of a mean from 2000 cells is 0.04 (lo) and 0.16 (hi)
    expect_equal(mean(r["lo", ]), 4, tolerance=0.05)
    expect_equal(mean(r["hi", ]), 50, tolerance=0.02)
    ## the rounding adds 1/12, relative error of the variance from 2000 draws is 3%
    expect_equal(var(r["lo", ]), 4, tolerance=0.2)
    expect_equal(var(r["hi", ]), 50, tolerance=0.2)
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend passes the table and dropout parameters of the object to the helper", {
    mat <- hs_const_matrix(ref_vals=c(1, 2, 4, 8, 16, 32), obs_vals=c(2, 3, 5, 9, 17, 33))
    obj <- hs_make_obj(mat, ref=list(n=1:5), obs=list(t=6:10))
    seen <- new.env()
    testthat::local_mocked_bindings(
        .get_simulated_cell_matrix_using_meanvar_trend_helper=function(gene_means, mean_var_table, num_cells, dropout_logistic_params=NULL) {
            seen$args <- list(gene_means, mean_var_table, num_cells, dropout_logistic_params)
            matrix(7, nrow=length(gene_means), ncol=num_cells)
        },
        .package="infercnv")
    gm <- c(a=1, b=2)
    r <- .get_simulated_cell_matrix_using_meanvar_trend(obj, gm, 3)
    expect_equal(r, matrix(7, 2, 3))
    expect_identical(seen$args[[1]], gm)
    expect_equal(seen$args[[2]], .get_mean_var_table(obj))
    expect_identical(seen$args[[3]], 3)
    expect_null(seen$args[[4]])

    hs_quiet(.get_simulated_cell_matrix_using_meanvar_trend(obj, gm, 3, include.dropout=TRUE))
    expect_equal(seen$args[[2]], .get_mean_var_table(obj))
    expect_true(is.list(seen$args[[4]]))
    expect_s3_class(seen$args[[4]]$spline, "smooth.spline")
    ## the dropout spline was fit to the observed (mean, p0) table of the object
    p0 <- hs_quiet(.get_logistic_params(.get_mean_vs_p0_table(obj)))
    expect_equal(predict(seen$args[[4]]$spline, log(c(2, 5, 20)))$y,
                 predict(p0$spline, log(c(2, 5, 20)))$y)
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend end to end reproduces gene means when there is no variance and no dropout", {
    mat <- hs_const_matrix(ref_vals=c(1, 2, 4, 8, 16, 32), obs_vals=c(2, 3, 5, 9, 17, 33))
    obj <- hs_make_obj(mat, ref=list(n=1:5), obs=list(t=6:10))
    gm <- c(g1=3.2, g2=5.1, g3=9.7, g4=1.4)
    expected <- matrix(c(3, 5, 10, 1), nrow=4, ncol=4,
                       dimnames=list(names(gm), paste0("sim_cell_", 1:4)))
    expect_equal(.get_simulated_cell_matrix_using_meanvar_trend(obj, gm, 4), expected)
    r <- hs_quiet(.get_simulated_cell_matrix_using_meanvar_trend(obj, gm, 4, include.dropout=TRUE))
    ## no zeros in the object, so the fitted dropout probability is 0 and nothing is dropped
    expect_equal(r, expected)
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix uses all cells unless groupings are given", {
    normal <- hs_const_matrix(ref_vals=c(1, 2, 4, 8, 16, 32), obs_vals=c(1, 2, 4, 8, 16, 32))
    gm <- c(g1=3.2, g2=5.1, g3=9.7, g4=1.4)
    expected <- matrix(c(3, 5, 10, 1), nrow=4, ncol=3,
                       dimnames=list(names(gm), paste0("sim_cell_", 1:3)))
    expect_equal(.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(
        gm, normal, 3, include.dropout=FALSE), expected)
    r <- hs_quiet(.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(
        gm, normal, 3, include.dropout=TRUE))
    expect_equal(r, expected)
    r <- hs_quiet(.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(
        gm, normal, 3, include.dropout=TRUE, cell_groupings=list(a=1:4, b=5:10)))
    expect_equal(r, expected)
})

test_that(".get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix hands the grouped tables to the helper", {
    normal <- matrix(c(1, 3, 0, 2, 6, 9,
                       2, 2, 4, 8, 0, 3,
                       5, 1, 1, 7, 2, 2,
                       0, 0, 3, 3, 3, 1), nrow=4, byrow=TRUE,
                     dimnames=list(paste0("g", 1:4), NULL))
    groups <- list(a=1:3, b=4:6)
    seen <- new.env()
    testthat::local_mocked_bindings(
        .get_simulated_cell_matrix_using_meanvar_trend_helper=function(gene_means, mean_var_table, num_cells, dropout_logistic_params=NULL) {
            seen$mvt <- mean_var_table
            seen$dropout <- dropout_logistic_params
            matrix(0, length(gene_means), num_cells)
        },
        .package="infercnv")
    .get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(
        c(a=1, b=2), normal, 2, include.dropout=FALSE, cell_groupings=groups)
    expect_equal(seen$mvt, .get_mean_var_given_matrix(normal, groups))
    expect_equal(seen$mvt$g, rep(c("a", "b"), each=4))
    expect_equal(seen$mvt$m[1:4], c(4/3, 8/3, 7/3, 1))
    expect_null(seen$dropout)
    hs_quiet(.get_simulated_cell_matrix_using_meanvar_trend_given_normal_matrix(
        c(a=1, b=2), normal, 2, include.dropout=TRUE, cell_groupings=groups))
    expect_s3_class(seen$dropout$spline, "smooth.spline")
})

## ---- simple simulation ----

test_that(".sim_expr_val returns 0 for a mean of 0 and draws from a negative binomial otherwise", {
    set.seed(1)
    before <- .Random.seed
    expect_identical(.sim_expr_val(0, NULL), 0)
    expect_identical(.Random.seed, before)
    for (m in c(0.5, 7, 80)) {
        set.seed(17)
        expected <- rnbinom(1, mu=m, size=1 / 0.25)
        set.seed(17)
        expect_identical(.sim_expr_val(m, NULL, common_dispersion=0.25), expected)
    }
    set.seed(2)
    d <- vapply(1:4000, function(i) .sim_expr_val(20, NULL, common_dispersion=0.1), 0)
    ## variance is mu + 0.1 mu^2 = 60, standard error of the mean 0.12
    expect_equal(mean(d), 20, tolerance=0.03)
    expect_equal(var(d), 60, tolerance=0.15)
})

test_that(".sim_expr_val dropout: spline and logistic forms", {
    sp1 <- hs_dropout_params(1)
    sp0 <- hs_dropout_params(0)
    expect_equal(.sim_expr_val(9, sp1), 0)
    expect_gt(.sim_expr_val(9, sp0), 0)
    ## with use_spline=FALSE the (midpt, slope) logistic in log(val) is used instead
    always_drop <- list(spline=hs_const_spline(0), midpt=-100, slope=1)
    never_drop <- list(spline=hs_const_spline(1), midpt=100, slope=1)
    set.seed(6)
    expect_equal(.sim_expr_val(9, always_drop, use_spline=FALSE), 0)
    expect_gt(.sim_expr_val(9, never_drop, use_spline=FALSE), 0)
    ## and the spline is what is used by default

})

test_that(".get_simulated_cell_matrix without a p0 table is a zero-aware negative binomial simulation", {
    gm <- c(a=20, b=0, c=3)
    set.seed(22)
    r <- .get_simulated_cell_matrix(gm, mean_p0_table=NULL, num_cells=1500, common_dispersion=0.1)
    expect_identical(dim(r), c(3L, 1500L))
    expect_identical(rownames(r), c("a", "b", "c"))
    expect_identical(colnames(r)[c(1, 1500)], c("spike_cell_1", "spike_cell_1500"))
    expect_true(all(r >= 0 & r == round(r)))
    expect_true(all(r["b", ] == 0))
    ## standard errors of the means: 0.2 and 0.05
    expect_equal(mean(r["a", ]), 20, tolerance=0.05)
    expect_equal(mean(r["c", ]), 3, tolerance=0.05)
})

test_that(".get_simulated_cell_matrix applies dropout derived from the p0 table", {
    ## p0 = 1 everywhere: the spline predicts 1 and every simulated value is dropped
    tb <- data.frame(m=c(1, 2, 4, 8, 16, 32, 64, 128), p0=1)
    gm <- c(a=20, b=5, c=3)
    r <- hs_quiet(.get_simulated_cell_matrix(gm, tb, num_cells=5, common_dispersion=0.1))
    expect_identical(dim(r), c(3L, 5L))
    expect_true(all(r == 0))
    ## p0 = 0 everywhere: nothing is dropped, so about all cells of a high-mean gene are positive
    tb0 <- data.frame(m=c(1, 2, 4, 8, 16, 32, 64, 128), p0=0)
    set.seed(9)
    r0 <- hs_quiet(.get_simulated_cell_matrix(c(a=200), tb0, num_cells=200, common_dispersion=0.1))
    expect_equal(mean(r0), 200, tolerance=0.1)
    expect_true(all(r0 > 0))
})

## ---- .build_and_add_hspike ----

hs_hspike_input <- function(with_ref=TRUE) {
    set.seed(1234)
    ng <- 30
    mat <- matrix(rnbinom(ng * 12, mu=rep(c(0, 0.5, 3, 10, 40), length.out=ng), size=2), nrow=ng,
                  dimnames=list(paste0("g", seq_len(ng)), paste0("c", 1:12)))
    ## guarantee zero-mean genes within each reference group
    mat[1, 1:8] <- 0
    mat[2, 5:8] <- 0
    if (with_ref) {
        hs_make_obj(mat, ref=list(A=1:4, B=5:8), obs=list(T1=9:12))
    } else {
        hs_make_obj(mat, ref=list(), obs=list(T1=9:12))
    }
}

## Replaces the matrix simulator with a deterministic one that records its calls.
hs_mock_meanvar <- function(env) {
    function(infercnv_obj, gene_means, num_cells, include.dropout=FALSE) {
        env$calls[[length(env$calls) + 1]] <- list(obj=infercnv_obj, gene_means=gene_means,
                                                   num_cells=num_cells, include.dropout=include.dropout)
        matrix(rep(round(gene_means * 1000), num_cells), nrow=length(gene_means), ncol=num_cells)
    }
}

test_that(".build_and_add_hspike builds a fake genome of eleven chromosomes and attaches the spike", {
    obj <- hs_hspike_input()
    env <- new.env(); env$calls <- list()
    testthat::local_mocked_bindings(.get_simulated_cell_matrix_using_meanvar_trend=hs_mock_meanvar(env),
                                    .package="infercnv")
    set.seed(11)
    out <- .build_and_add_hspike(obj, "meanvar")
    hs <- out@.hspike
    expect_s4_class(hs, "infercnv")
    ## the host object is otherwise unchanged
    expect_identical(out@expr.data, obj@expr.data)
    expect_identical(out@reference_grouped_cell_indices, obj@reference_grouped_cell_indices)

    ## 30 input genes: 10 * 400 + max(30 - 4000, 400) = 4400 fake genes
    expect_identical(dim(hs@expr.data), c(4400L, 400L))
    expect_equal(rownames(hs@expr.data), paste0("gene_", 1:4400))
    expect_equal(rownames(hs@gene_order), paste0("gene_", 1:4400))
    go <- hs@gene_order
    expect_equal(as.character(go$chr),
                 rep(names(hs_cnv_by_chr), times=c(rep(400, 10), 400)))
    expect_equal(go$start, rep(1:400, 11))
    expect_equal(go$stop, go$start)

    ## per group: 100 normal cells then 100 spiked tumour cells
    expect_equal(names(hs@reference_grouped_cell_indices), c("simnorm_cell_A", "simnorm_cell_B"))
    expect_equal(names(hs@observation_grouped_cell_indices), c("spike_tumor_cell_A", "spike_tumor_cell_B"))
    expect_equal(hs@reference_grouped_cell_indices$simnorm_cell_A, 1:100)
    expect_equal(hs@observation_grouped_cell_indices$spike_tumor_cell_A, 101:200)
    expect_equal(hs@reference_grouped_cell_indices$simnorm_cell_B, 201:300)
    expect_equal(hs@observation_grouped_cell_indices$spike_tumor_cell_B, 301:400)
    expect_equal(colnames(hs@expr.data)[c(1, 100, 101, 200, 201, 301, 400)],
                 c("simnorm_cell_A1", "simnorm_cell_A100", "spike_tumor_cell_A1", "spike_tumor_cell_A100",
                   "simnorm_cell_B1", "spike_tumor_cell_B1", "spike_tumor_cell_B100"))
    expect_identical(hs@count.data, hs@expr.data * 0 + hs@count.data)
})

test_that(".build_and_add_hspike draws gene means from the normal cells and multiplies them by the chromosome copy number", {
    obj <- hs_hspike_input()
    env <- new.env(); env$calls <- list()
    testthat::local_mocked_bindings(.get_simulated_cell_matrix_using_meanvar_trend=hs_mock_meanvar(env),
                                    .package="infercnv")
    set.seed(11)
    out <- .build_and_add_hspike(obj, "meanvar")
    ## two reference groups, normals then spiked each: four simulator calls
    expect_length(env$calls, 4)
    for (cl in env$calls) {
        expect_identical(cl$num_cells, 100)
        expect_true(cl$include.dropout)
        expect_equal(names(cl$gene_means), paste0("gene_", 1:4400))
        expect_s4_class(cl$obj, "infercnv")
    }

    ## expected gene indices: the first random draw after the seed
    set.seed(11)
    idx <- sample(seq_len(30), size=4400, replace=TRUE)
    cnv_per_gene <- rep(unname(hs_cnv_by_chr), each=400)
    for (grp in list(list(call=1, cols=1:4), list(call=3, cols=5:8))) {
        orig <- rowMeans(obj@expr.data[, grp$cols])[idx]
        orig[orig == 0] <- 1e-3
        expect_equal(unname(env$calls[[grp$call]]$gene_means), unname(orig))
        expect_equal(unname(env$calls[[grp$call + 1]]$gene_means), unname(orig * cnv_per_gene))
    }
    ## zero-mean genes were replaced by 1e-3, and the draw includes at least one such gene
    expect_true(any(env$calls[[1]]$gene_means == 1e-3))
    ## the hspike matrix is normalised so each cell sums to the median normal library size of the last group
    target <- median(colSums(obj@expr.data[, 5:8]))
    expect_equal(unname(colSums(out@.hspike@expr.data)), rep(target, 400))
    ## raw simulated values are gene_mean * 1000: proportions are preserved by the normalisation
    raw1 <- round(env$calls[[1]]$gene_means * 1000)
    expect_equal(unname(out@.hspike@expr.data[, 1]), unname(raw1 / sum(raw1) * target))
    raw2 <- round(env$calls[[2]]$gene_means * 1000)
    expect_equal(unname(out@.hspike@expr.data[, 101]), unname(raw2 / sum(raw2) * target))
})

test_that(".build_and_add_hspike with aggregate_normals pools the reference groups into one", {
    obj <- hs_hspike_input()
    env <- new.env(); env$calls <- list()
    testthat::local_mocked_bindings(.get_simulated_cell_matrix_using_meanvar_trend=hs_mock_meanvar(env),
                                    .package="infercnv")
    set.seed(11)
    out <- .build_and_add_hspike(obj, "meanvar", aggregate_normals=TRUE)
    hs <- out@.hspike
    expect_length(env$calls, 2)
    expect_identical(dim(hs@expr.data), c(4400L, 200L))
    expect_equal(names(hs@reference_grouped_cell_indices), "simnorm_cell_normalsToUse")
    expect_equal(names(hs@observation_grouped_cell_indices), "spike_tumor_cell_normalsToUse")
    expect_equal(colnames(hs@expr.data)[c(1, 101)], c("simnorm_cell_normalsToUse1", "spike_tumor_cell_normalsToUse1"))
    set.seed(11)
    idx <- sample(seq_len(30), size=4400, replace=TRUE)
    orig <- rowMeans(obj@expr.data[, 1:8])[idx]
    orig[orig == 0] <- 1e-3
    expect_equal(unname(env$calls[[1]]$gene_means), unname(orig))
    expect_equal(unname(colSums(hs@expr.data)), rep(median(colSums(obj@expr.data[, 1:8])), 200))
})

test_that(".build_and_add_hspike without reference cells uses the observation cells as the normals", {
    obj <- hs_hspike_input(with_ref=FALSE)
    env <- new.env(); env$calls <- list()
    testthat::local_mocked_bindings(.get_simulated_cell_matrix_using_meanvar_trend=hs_mock_meanvar(env),
                                    .package="infercnv")
    set.seed(11)
    out <- .build_and_add_hspike(obj, "meanvar")
    hs <- out@.hspike
    expect_length(env$calls, 2)
    expect_equal(names(hs@reference_grouped_cell_indices), "simnorm_cell_normalsToUse")
    expect_equal(names(hs@observation_grouped_cell_indices), "spike_tumor_cell_normalsToUse")
    expect_identical(dim(hs@expr.data), c(4400L, 200L))
    ## the simulator sees an object whose observation grouping is the proxy-normal grouping
    for (cl in env$calls) {
        expect_equal(names(cl$obj@observation_grouped_cell_indices), "normalsToUse")
        expect_equal(unname(cl$obj@observation_grouped_cell_indices$normalsToUse), 9:12)
    }
    set.seed(11)
    idx <- sample(seq_len(30), size=4400, replace=TRUE)
    orig <- rowMeans(obj@expr.data[, 9:12])[idx]
    orig[orig == 0] <- 1e-3
    expect_equal(unname(env$calls[[1]]$gene_means), unname(orig))
    ## the returned object keeps its original grouping
    expect_equal(names(out@observation_grouped_cell_indices), "T1")
    expect_equal(length(out@reference_grouped_cell_indices), 0)
    expect_equal(unname(colSums(hs@expr.data)), rep(median(colSums(obj@expr.data[, 9:12])), 200))
})

test_that(".build_and_add_hspike gives the last fake chromosome the remaining genes of large inputs", {
    set.seed(5)
    ng <- 4700
    mat <- matrix(rpois(ng * 6, 4), nrow=ng, dimnames=list(paste0("g", seq_len(ng)), paste0("c", 1:6)))
    obj <- hs_make_obj(mat, ref=list(R=1:3), obs=list(T=4:6))
    env <- new.env(); env$calls <- list()
    testthat::local_mocked_bindings(.get_simulated_cell_matrix_using_meanvar_trend=hs_mock_meanvar(env),
                                    .package="infercnv")
    out <- .build_and_add_hspike(obj, "meanvar")
    hs <- out@.hspike
    ## 4700 - 4000 = 700 remaining genes on chrF
    expect_identical(dim(hs@expr.data), c(4700L, 200L))
    expect_equal(sum(hs@gene_order$chr == "chr_F"), 700)
    expect_equal(sum(hs@gene_order$chr == "chr_3pt0"), 400)
    expect_equal(hs@gene_order$start[hs@gene_order$chr == "chr_F"], 1:700)
})

test_that(".build_and_add_hspike meanvar run end to end gives normal cells and copy-number-scaled spiked cells", {
    set.seed(3)
    mu <- rep(c(0.5, 3, 10, 40), 15)
    mat <- matrix(rnbinom(60 * 20, mu=mu, size=2), 60, 20,
                  dimnames=list(paste0("g", 1:60), paste0("c", 1:20)))
    obj <- hs_make_obj(mat, ref=list(n=1:10), obs=list(t=11:20))
    set.seed(41)
    out <- hs_quiet(.build_and_add_hspike(obj, "meanvar"))
    hs <- out@.hspike
    expect_identical(dim(hs@expr.data), c(4400L, 200L))
    expect_true(all(hs@expr.data >= 0))
    target <- median(colSums(mat[, 1:10]))
    ## every cell with any counts is scaled to the reference library size
    expect_equal(unname(colSums(hs@expr.data)), rep(target, 200))
    chr <- hs@gene_order$chr
    spiked <- hs@expr.data[, hs@observation_grouped_cell_indices[[1]]]
    normal <- hs@expr.data[, hs@reference_grouped_cell_indices[[1]]]
    chr_means <- function(m) tapply(rowMeans(m), factor(chr, levels=names(hs_cnv_by_chr)), mean)
    sm <- chr_means(spiked)
    nm <- chr_means(normal)
    ## a 0.01-fold region is almost empty compared with the neutral regions; the 3-fold region is high
    expect_lt(sm[["chr_0"]], 0.25 * mean(sm[c("chrA", "chr_B", "chr_C")]))
    expect_gt(sm[["chr_3pt0"]], 1.5 * mean(sm[c("chr_D", "chr_E", "chr_F")]))
    expect_lt(sm[["chr_0pt5"]], sm[["chr_1pt5"]])
    ## normal cells show no such structure
    expect_lt(max(nm) / min(nm), 1.5)
})

test_that(".build_and_add_hspike simple run end to end returns a normalised spike with the same layout", {
    set.seed(3)
    mu <- rep(c(0.5, 3, 10, 40), 15)
    mat <- matrix(rnbinom(60 * 20, mu=mu, size=2), 60, 20,
                  dimnames=list(paste0("g", 1:60), paste0("c", 1:20)))
    obj <- hs_make_obj(mat, ref=list(n=1:10), obs=list(t=11:20))
    set.seed(42)
    out <- hs_quiet(.build_and_add_hspike(obj, "simple"))
    hs <- out@.hspike
    expect_identical(dim(hs@expr.data), c(4400L, 200L))
    expect_equal(names(hs@reference_grouped_cell_indices), "simnorm_cell_n")
    expect_equal(names(hs@observation_grouped_cell_indices), "spike_tumor_cell_n")
    expect_true(all(hs@expr.data >= 0))
    expect_equal(unname(colSums(hs@expr.data)), rep(median(colSums(mat[, 1:10])), 200))
    chr <- hs@gene_order$chr
    spiked <- rowMeans(hs@expr.data[, 101:200])
    expect_lt(mean(spiked[chr == "chr_0"]), 0.25 * mean(spiked[chr == "chrA"]))
    expect_gt(mean(spiked[chr == "chr_3pt0"]), 1.5 * mean(spiked[chr == "chr_E"]))
})

test_that(".build_and_add_hspike rejects an unknown simulation method", {
    expect_error(.build_and_add_hspike(hs_hspike_input(), "nonsense"), "should be one of")
})

## ---- .sim_foreground ----

test_that(".sim_foreground replaces every cell by simulated values and rescales to the reference library size (meanvar)", {
    ref_vals <- c(2, 4, 6, 10, 20, 30)
    obs_vals <- c(5, 1, 9, 12, 3, 40)
    mat <- hs_const_matrix(ref_vals, obs_vals)
    obj <- hs_make_obj(mat, ref=list(n=1:5), obs=list(t=6:10))
    out <- hs_quiet(.sim_foreground(obj, "meanvar"))
    ## constant genes: no variance and no dropout, so the simulated cell is the (rounded) group mean
    target <- median(colSums(mat[, 1:5]))
    expect_equal(target, sum(ref_vals))
    expected_ref <- ref_vals / sum(ref_vals) * target
    expected_obs <- obs_vals / sum(obs_vals) * target
    expect_equal(unname(out@expr.data[, 1]), expected_ref)
    expect_equal(unname(out@expr.data[, 5]), expected_ref)
    expect_equal(unname(out@expr.data[, 6]), expected_obs)
    expect_equal(unname(out@expr.data[, 10]), expected_obs)
    expect_identical(dimnames(out@expr.data), dimnames(mat))
    ## everything but expr.data is untouched
    expect_identical(out@count.data, obj@count.data)
    expect_identical(out@gene_order, obj@gene_order)
})

test_that(".sim_foreground simulates each sample with its own gene means", {
    env <- new.env(); env$calls <- list()
    mat <- matrix(c(1, 3, 10, 20,
                    2, 2,  0,  0,
                    4, 6, 30, 50), nrow=3, byrow=TRUE,
                  dimnames=list(c("a", "b", "c"), paste0("c", 1:4)))
    obj <- hs_make_obj(mat, ref=list(R=1:2), obs=list(T=3:4))
    testthat::local_mocked_bindings(
        .get_simulated_cell_matrix_using_meanvar_trend=function(infercnv_obj, gene_means, num_cells, include.dropout=FALSE) {
            env$calls[[length(env$calls) + 1]] <- list(gene_means=gene_means, num_cells=num_cells, drop=include.dropout)
            matrix(seq_along(gene_means), length(gene_means), num_cells)
        }, .package="infercnv")
    out <- .sim_foreground(obj, "meanvar")
    ## samples are visited observations first, then references
    expect_length(env$calls, 2)
    expect_equal(env$calls[[1]]$gene_means, c(a=15, b=1e-3, c=40))
    expect_equal(env$calls[[2]]$gene_means, c(a=2, b=2, c=5))
    expect_true(all(vapply(env$calls, function(x) x$drop, NA)))
    expect_true(all(vapply(env$calls, function(x) x$num_cells, 0) == 2))
    ## values 1,2,3 per cell scaled to the target library size (median colsum of reference = 9)
    target <- median(colSums(mat[, 1:2]))
    expect_equal(unname(out@expr.data[, 1]), c(1, 2, 3) / 6 * target)
    expect_equal(unname(out@expr.data[, 4]), c(1, 2, 3) / 6 * target)
})

test_that(".sim_foreground with the simple method keeps shape and library size", {
    set.seed(8)
    mat <- matrix(rnbinom(20 * 12, mu=rep(c(8, 25, 60, 120), 5), size=5), 20, 12,
                  dimnames=list(paste0("g", 1:20), paste0("c", 1:12)))
    mat[mat == 0] <- 1  # no zeros: the fitted dropout is 0
    obj <- hs_make_obj(mat, ref=list(n=1:6), obs=list(t=7:12))
    set.seed(9)
    out <- hs_quiet(.sim_foreground(obj, "simple"))
    expect_identical(dim(out@expr.data), dim(mat))
    expect_identical(dimnames(out@expr.data), dimnames(mat))
    expect_true(all(out@expr.data >= 0))
    expect_equal(unname(colSums(out@expr.data)), rep(median(colSums(mat[, 1:6])), 12))
    expect_false(isTRUE(all.equal(out@expr.data, mat)))
    ## the simulated data keeps the gene ranking of the originals on average
    expect_gt(cor(rowMeans(out@expr.data), rowMeans(mat)), 0.9)
})

## ---- SplatterScrape ----

test_that(".winsorize clamps values to the q and 1-q quantiles", {
    x <- c(1, 2, 3, 4, 100, -50, 5, 6, 7, 8)
    lohi <- quantile(x, c(0.1, 0.9))
    expect_equal(lohi[[1]], -4.1 * 1 + 0, tolerance=1e-12)
    expect_equal(.winsorize(x, 0.1), c(1, 2, 3, 4, 17.2, -4.1, 5, 6, 7, 8))
    ## missing values are left in place
    expect_equal(.winsorize(c(1, NA, 2, 3, 4, 5, 6, 7, 8, 9, 10), 0.1)[2], NA_real_)
    ## q = 0 leaves the data alone
    expect_equal(.winsorize(x, 0), x)
    ## q > 0.5 reverses the quantile pair
    expect_equal(.winsorize(c(0, 5, 10), 0.75), c(2.5, 5, 7.5))
})

test_that(".logistic is the logistic function", {
    expect_equal(.logistic(1, x0=1, k=5), 0.5)
    expect_equal(.logistic(c(0, 2), x0=1, k=1), c(1 / (1 + exp(1)), 1 / (1 + exp(-1))))
    expect_equal(.logistic(3, x0=1, k=-2), 1 / (1 + exp(4)))
})

hs_splat_params <- function(nGenes=20, nCells=30, ...) {
    p <- list(include.dropout=FALSE, use.spline.dropout.fit=FALSE,
              mean.shape=2, mean.rate=0.5,
              lib.loc=8, lib.scale=0.2, lib.norm=FALSE,
              out.prob=0, out.facLoc=1, out.facScale=0.5,
              bcv.common=0.1, bcv.df=60,
              dropout.mid=0, dropout.shape=-1,
              nGenes=nGenes, nCells=nCells)
    extra <- list(...)
    p[names(extra)] <- extra
    p
}

test_that(".getLNormFactors selects no factors with probability 0 and only factors >= 1 when never reversed", {
    set.seed(1)
    expect_equal(.getLNormFactors(50, 0, 0, 1, 0.5), rep(1, 50))
    set.seed(2)
    f <- .getLNormFactors(500, 0.4, 0, 1, 0.5)
    expect_length(f, 500)
    ## neg.prob = 0: factors < 1 get their direction flipped, so every factor is >= 1
    expect_true(all(f >= 1))
    ## selection fraction 0.4 of 500: sd 0.022
    expect_equal(mean(f != 1), 0.4, tolerance=0.1)
    set.seed(3)
    f_all <- .getLNormFactors(500, 1, 0, 0, 1)
    expect_true(all(f_all >= 1))
    expect_gt(mean(f_all > 1), 0.99)
    ## neg.prob = 1: directions reversed, so all factors <= 1
    set.seed(4)
    expect_true(all(.getLNormFactors(500, 1, 1, 0, 1) <= 1))
})

test_that(".splatEstLib fits a normal to normally distributed library sizes and a log-normal otherwise", {
    skip_if_not_installed("fitdistrplus")
    set.seed(5)
    cnt_norm <- matrix(0, 3, 100)
    lib <- rnorm(100, 1e4, 500)
    cnt_norm[1, ] <- lib
    expect_gt(shapiro.test(lib)$p.value, 0.2)
    r <- suppressMessages(.splatEstLib(cnt_norm, list()))
    expect_true(r$lib.norm)
    expect_equal(r$lib.loc, mean(lib), tolerance=1e-4)
    expect_equal(r$lib.scale, sqrt(mean((lib - mean(lib))^2)), tolerance=1e-3)

    set.seed(6)
    lib2 <- exp(rnorm(60, 9, 1))
    expect_lt(shapiro.test(lib2)$p.value, 0.2)
    cnt_ln <- matrix(0, 3, 60)
    cnt_ln[1, ] <- lib2
    r2 <- .splatEstLib(cnt_ln, list(keep="me"))
    expect_false(r2$lib.norm)
    expect_equal(r2$keep, "me")
    expect_equal(r2$lib.loc, mean(log(lib2)), tolerance=1e-4)
    expect_equal(r2$lib.scale, sqrt(mean((log(lib2) - mean(log(lib2)))^2)), tolerance=1e-3)
})

test_that(".splatEstMean fits a gamma to the winsorised nonzero means", {
    skip_if_not_installed("fitdistrplus")
    set.seed(7)
    means <- rgamma(400, shape=3, rate=1.5)
    norm.counts <- matrix(rep(means, 5), ncol=5)
    r <- suppressWarnings(.splatEstMean(norm.counts, list(a=1)))
    expect_equal(r$a, 1)
    expect_gt(r$mean.shape, 0)
    expect_gt(r$mean.rate, 0)
    ## fitted gamma mean is the mean of the winsorised data to within a few percent
    wm <- mean(.winsorize(means, 0.1))
    expect_equal(r$mean.shape / r$mean.rate, wm, tolerance=0.1)
})

test_that(".splatEstOutlier counts the share of genes with log-mean above median + 2 MAD", {
    skip_if_not_installed("fitdistrplus")
    set.seed(8)
    lm <- c(rnorm(95, 1, 0.3), rep(4, 5))
    means <- exp(lm)
    norm.counts <- matrix(rep(means, 3), ncol=3)
    r <- .splatEstOutlier(norm.counts, list())
    expected_bound <- median(lm) + 2 * mad(lm)
    outs <- which(lm > expected_bound)
    expect_equal(r$out.prob, length(outs) / 100)
    expect_gte(length(outs), 5)
    facs <- means[outs] / median(means)
    expect_equal(r$out.facLoc, mean(log(facs)), tolerance=1e-3)
    expect_equal(r$out.facScale, sqrt(mean((log(facs) - mean(log(facs)))^2)), tolerance=1e-3)
    ## without outliers only the probability is set
    r0 <- .splatEstOutlier(matrix(rep(exp(seq(0, 1, length.out=100)), 3), ncol=3), list())
    expect_equal(r0$out.prob, 0)
    expect_null(r0$out.facLoc)
})

test_that(".splatEstBCV maps edgeR's common dispersion to bcv.common", {
    skip_if_not_installed("edgeR")
    set.seed(9)
    cnt <- matrix(rnbinom(200 * 30, mu=rep(c(10, 30, 80, 150), 50), size=2), nrow=200)
    r <- .splatEstBCV(cnt, list(z=0))
    d <- edgeR::estimateDisp(cnt, design=matrix(1, 30, 1))
    expect_equal(r$bcv.common, 0.1 + 0.25 * d$common.dispersion)
    expect_equal(r$bcv.df, d$prior.df)
    expect_equal(r$z, 0)
})

test_that(".splatEstDropout fits the logistic and spline to the zero fraction against log mean", {
    set.seed(10)
    x <- seq(-2, 4, length.out=80)
    ncell <- 100
    ## per-gene zero fraction follows a logistic in log mean
    pz <- 1 / (1 + exp(1.2 * (x - 1)))
    norm.counts <- t(vapply(seq_along(x), function(i) {
        v <- rep(exp(x[i]), ncell)
        v[seq_len(min(round(pz[i] * ncell), ncell - 1))] <- 0
        v * exp(x[i]) * ncell / sum(v)  # rescale so the row mean is exactly exp(x[i])
    }, numeric(ncell)))
    means <- rowMeans(norm.counts)
    r <- suppressWarnings(.splatEstDropout(norm.counts, list(keep=TRUE)))
    expect_true(r$keep)
    y <- rowSums(norm.counts == 0) / ncell
    xs <- log(means)
    sp <- predict(r$dropout.spline.fit, xs)$y
    expect_equal(sp, y, tolerance=0.1)
    expect_lt(r$dropout.shape, 0)
    expect_equal(r$dropout.shape, -1.2, tolerance=0.3)
    expect_equal(r$dropout.mid, 1, tolerance=0.5)
    ## the logistic fit reproduces the zero fractions
    expect_equal(.logistic(xs, r$dropout.mid, r$dropout.shape), y, tolerance=0.1)
})

test_that(".simulateSingleCellCountsMatrixSplatterScrape returns counts with Splat's names and shape", {
    skip_if_not_installed("SingleCellExperiment")
    p <- hs_splat_params(nGenes=20, nCells=30)
    set.seed(11)
    sim <- suppressMessages(.simulateSingleCellCountsMatrixSplatterScrape(p))
    cnt <- BiocGenerics::counts(sim)
    expect_identical(dim(cnt), c(20L, 30L))
    expect_equal(rownames(cnt), paste0("Gene", 1:20))
    expect_equal(colnames(cnt), paste0("Cell", 1:30))
    expect_true(all(cnt >= 0 & cnt == round(cnt)))
    ## without dropout the counts equal the true counts
    expect_equal(SummarizedExperiment::assay(sim, "TrueCounts"), cnt)
})

test_that(".simulateSingleCellCountsMatrixSplatterScrape rejects use.genes.means of the wrong length", {
    skip_if_not_installed("SingleCellExperiment")
    p <- hs_splat_params(nGenes=20, nCells=5)
    expect_error(.simulateSingleCellCountsMatrixSplatterScrape(p, use.genes.means=rep(1, 19)),
                 "not matching the params nGenes")
})

test_that("splat simulation of given gene means keeps their ranking and library sizes", {
    skip_if_not_installed("SingleCellExperiment")
    ## gene means are used as is (no outliers); expected library size is exp(lib.loc + lib.scale^2/2)
    gm <- rep(c(1, 10, 100, 1000), each=10)
    p <- hs_splat_params(nGenes=40, nCells=300, lib.loc=9, lib.scale=0.1)
    set.seed(12)
    sim <- suppressMessages(.simulateSingleCellCountsMatrixSplatterScrape(p, use.genes.means=gm))
    cnt <- BiocGenerics::counts(sim)
    expect_equal(unname(SummarizedExperiment::rowData(sim)$GeneMean), gm)
    expect_equal(unname(SummarizedExperiment::rowData(sim)$BaseGeneMean), gm)
    expect_equal(unname(SummarizedExperiment::rowData(sim)$OutlierFactor), rep(1, 40))
    gmean <- rowMeans(cnt)
    ## group means follow the requested relative abundance: group 4 has 10x the mean of group 3
    gmeans <- tapply(gmean, rep(1:4, each=10), mean)
    expect_equal(unname(gmeans[2] / gmeans[1]), 10, tolerance=0.4)
    expect_equal(unname(gmeans[4] / gmeans[3]), 10, tolerance=0.3)
    ## mean library size close to exp(9 + 0.005) = 8140
    expect_equal(mean(colSums(cnt)), exp(9 + 0.1^2 / 2), tolerance=0.1)
})

test_that("splat dropout probability follows the logistic in log(cell mean)", {
    skip_if_not_installed("SingleCellExperiment")
    gm <- rep(c(1, 10, 100, 1000), each=10)
    p <- hs_splat_params(nGenes=40, nCells=50, include.dropout=TRUE, dropout.mid=3, dropout.shape=-1.5)
    set.seed(13)
    sim <- suppressMessages(utils::capture.output(
        s <- .simulateSingleCellCountsMatrixSplatterScrape(p, use.genes.means=gm)))
    dp <- SummarizedExperiment::assay(s, "DropProb")
    cm <- SummarizedExperiment::assay(s, "CellMeans")
    expect_equal(unname(dp), unname(1 / (1 + exp(1.5 * (log(cm) - 3)))), tolerance=1e-8)
    expect_true(all(dp >= 0 & dp <= 1))
    cnt <- BiocGenerics::counts(s)
    dropped <- SummarizedExperiment::assay(s, "Dropout")
    expect_true(all(cnt[dropped] == 0))
    tc <- SummarizedExperiment::assay(s, "TrueCounts")
    expect_equal(cnt[!dropped], tc[!dropped])
})

test_that("splat dropout with a spline fit clamps probabilities to [0, 1]", {
    skip_if_not_installed("SingleCellExperiment")
    gm <- rep(c(5, 50, 500, 5000), each=5)
    ## a spline that is 0.3 everywhere, so the dropout probability is 0.3 for any mean
    sp <- hs_const_spline(0.3)
    sp$fit$coef <- sp$fit$coef
    p <- hs_splat_params(nGenes=20, nCells=400, include.dropout=TRUE, use.spline.dropout.fit=TRUE,
                         dropout.spline.fit=sp)
    set.seed(14)
    utils::capture.output(s <- suppressMessages(
        .simulateSingleCellCountsMatrixSplatterScrape(p, use.genes.means=gm)))
    dp <- SummarizedExperiment::assay(s, "DropProb")
    expect_equal(unname(dp), matrix(0.3, 20, 400), tolerance=1e-6)
    ## standard deviation of the dropped fraction over 8000 draws is 0.005
    expect_equal(mean(SummarizedExperiment::assay(s, "Dropout")), 0.3, tolerance=0.05)
})

test_that(".estimateSingleCellParamsSplatterScrape output feeds the splat simulator for the h-spike path", {
    skip_if_not_installed("SingleCellExperiment")
    skip_if_not_installed("fitdistrplus")
    skip_if_not_installed("edgeR")
    set.seed(15)
    mu <- rgamma(300, shape=2, rate=0.5)
    cnt <- matrix(rnbinom(300 * 50, mu=mu, size=5), 300, 50,
                  dimnames=list(paste0("g", 1:300), paste0("c", 1:50)))
    pr <- NULL
    utils::capture.output(pr <- suppressWarnings(suppressMessages(.estimateSingleCellParamsSplatterScrape(cnt))))
    expect_identical(pr$nGenes, 300L)
    expect_identical(pr$nCells, 50L)
    expect_false(pr$include.dropout)
    expect_false(pr$use.spline.dropout.fit)
    expect_true(all(c("mean.shape", "mean.rate", "lib.loc", "lib.scale", "lib.norm", "out.prob",
                      "bcv.common", "bcv.df", "dropout.mid", "dropout.shape", "dropout.spline.fit") %in% names(pr)))
    expect_gt(pr$mean.shape, 0)
    expect_gt(pr$mean.rate, 0)
    ## library size of the simulated cells is on the scale of the input
    expect_equal(pr$lib.loc, mean(colSums(cnt)), tolerance=0.05)
    pr$nGenes <- 100
    pr$nCells <- 20
    pr$out.facLoc <- 1
    pr$out.facScale <- 0.5
    set.seed(16)
    s <- suppressWarnings(suppressMessages(.simulateSingleCellCountsMatrixSplatterScrape(pr, use.genes.means=rep(1:5, 20))))
    expect_identical(dim(BiocGenerics::counts(s)), c(100L, 20L))
})
