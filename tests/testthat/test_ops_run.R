# Tests for the per-step operations that run() calls (R/inferCNV_ops.R), noise reduction
# (R/noise_reduction.R), tumor subclustering (R/inferCNV_tumor_subclusters*.R), non-DE masking
# (R/inferCNV_mask_non_DE.R) and a tiny end-to-end run().

futile.logger::flog.threshold(futile.logger::ERROR)
assign("GLOBAL_NUM_THREADS", 1, envir=infercnv:::infercnv.env)

## ---------------------------------------------------------------------------------------------
## helpers (prefix opsrun_ so that sibling test files cannot collide)
## ---------------------------------------------------------------------------------------------

## genes in rows, cells in columns; chr is one entry per gene
opsrun_obj <- function(expr, chr=rep("chr1", nrow(expr)), ref=NULL, obs=NULL, start=NULL) {
	if (is.null(rownames(expr))) rownames(expr) <- paste0("g", seq_len(nrow(expr)))
	if (is.null(colnames(expr))) colnames(expr) <- paste0("c", seq_len(ncol(expr)))
	if (is.null(start)) start <- seq_len(nrow(expr)) * 1000
	gene_order <- data.frame(chr=factor(chr, levels=unique(chr)), start=start, stop=start + 500,
							 row.names=rownames(expr))
	if (is.null(obs)) obs <- list(tumor=setdiff(seq_len(ncol(expr)), unlist(ref)))
	if (is.null(ref)) ref <- list()
	methods::new("infercnv", expr.data=expr, count.data=expr, gene_order=gene_order,
				 reference_grouped_cell_indices=ref, observation_grouped_cell_indices=obs,
				 tumor_subclusters=NULL, options=list(), .hspike=NULL)
}

## reference implementation of the pyramid-weighted moving average with renormalised, truncated
## windows at the ends; NAs are removed first and put back afterwards
opsrun_ref_pyramid <- function(x, w) {
	nas <- is.na(x)
	v <- x[!nas]
	n <- length(v)
	tl <- (w - 1) / 2
	out <- numeric(n)
	for (i in seq_len(n)) {
		k <- max(1, i - tl):min(n, i + tl)
		wt <- tl + 1 - abs(k - i)
		out[i] <- sum(v[k] * wt) / sum(wt)
	}
	x[!nas] <- out
	x
}

## reference for caTools::runmean(k, endrule="mean"): equal weights, window truncated at the ends
opsrun_ref_runmean <- function(x, k) {
	h <- (k - 1) / 2
	n <- length(x)
	vapply(seq_len(n), function(i) mean(x[max(1, i - h):min(n, i + h)]), numeric(1))
}

opsrun_with_hspike <- function(obj, hspike) {
	obj@.hspike <- hspike
	obj
}

## ---------------------------------------------------------------------------------------------
## smoothing
## ---------------------------------------------------------------------------------------------

test_that(".smooth_center_helper applies the pyramid filter (1,2,1)/4 and keeps the ends", {
	expect_equal(infercnv:::.smooth_center_helper(c(1, 2, 3, 4, 5, 6), 3),
				 c(1, 2, 3, 4, 5, 6))
	## hand: position 2 = (1*1 + 2*2 + 1*4)/4 = 9/4, position 3 = (2 + 8 + 8)/4, position 4 = (4 + 16 + 16)/4
	expect_equal(infercnv:::.smooth_center_helper(c(1, 2, 4, 8, 16), 3),
				 c(1, 9/4, 4.5, 9, 16))
	## width 5, weights (1,2,3,2,1)/9
	expect_equal(infercnv:::.smooth_center_helper(c(1, 2, 4, 8, 16, 32, 64), 5),
				 c(1, 2, 49/9, 98/9, 196/9, 32, 64))
})

test_that(".smooth_center_helper skips NA positions and smooths across the remaining values", {
	## non-NA values 1,4,8,16 -> interior (1+8+8)/4 = 4.25 and (4+16+16)/4 = 9
	expect_equal(infercnv:::.smooth_center_helper(c(1, NA, 4, 8, 16), 3),
				 c(1, NA, 4.25, 9, 16))
})

test_that(".smooth_helper uses renormalised truncated windows at the ends (hand-worked)", {
	## window 3, x = 1,2,4,8,16
	## end 1: weights (2,1) over (1,2) = (2+2)/3 ; end 5: weights (1,2) over (8,16) = (8+32)/3
	expect_equal(infercnv:::.smooth_helper(c(1, 2, 4, 8, 16), 3),
				 c(4/3, 9/4, 4.5, 9, 40/3))
	## window 5 longer than data (3 values): weights at offsets 0,1,2 = 3,2,1
	## pos 1 = (1*3 + 2*2 + 4*1)/6 ; pos 2 = (1*2 + 2*3 + 4*2)/7 ; pos 3 = (4*3 + 2*2 + 1*1)/6
	expect_equal(infercnv:::.smooth_helper(c(1, 2, 4), 5),
				 c(11/6, 16/7, 17/6))
	## a single value is returned unchanged
	expect_equal(infercnv:::.smooth_helper(7, 3), 7)
})

test_that(".smooth_helper matches the reference pyramid average for many lengths and windows", {
	for (n in c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 20)) {
		for (w in c(3, 5, 7, 11)) {
			x <- sin(seq_len(n)) * 10 + seq_len(n) %% 3
			expect_equal(infercnv:::.smooth_helper(x, w), opsrun_ref_pyramid(x, w),
						 info=sprintf("n=%d window=%d", n, w))
		}
	}
})

test_that(".smooth_helper keeps NA positions and smooths the remaining values", {
	x <- c(1, NA, 4, 8, 16, NA, 3, 9, 2, 7)
	res <- infercnv:::.smooth_helper(x, 3)
	expect_equal(is.na(res), is.na(x))
	expect_equal(res, opsrun_ref_pyramid(x, 3))
	## hand: first non-NA value uses weights (2,1) over (1,4): (2 + 4)/3 = 2
	expect_equal(res[1], 2)
	## hand: third entry (second non-NA, value 4) = (1*1 + 2*4 + 1*8)/4 = 4.25
	expect_equal(res[3], 4.25)
})

test_that(".smooth_window smooths columns independently and keeps dimnames", {
	m <- cbind(a=c(1, 2, 4, 8, 16), b=c(16, 8, 4, 2, 1))
	rownames(m) <- paste0("g", 1:5)
	res <- infercnv:::.smooth_window(m, 3)
	expect_equal(dimnames(res), dimnames(m))
	expect_equal(res[, "a"], setNames(c(4/3, 9/4, 4.5, 9, 40/3), rownames(m)))
	expect_equal(res[, "b"], setNames(rev(c(4/3, 9/4, 4.5, 9, 40/3)), rownames(m)))
})

test_that(".smooth_window returns the data unchanged for window_length < 2", {
	m <- matrix(c(3, 1, 4, 1, 5, 9), ncol=2)
	expect_identical(infercnv:::.smooth_window(m, 1), m)
	expect_identical(infercnv:::.smooth_window(m, 0), m)
})

test_that("smooth_by_chromosome smooths each chromosome separately and keeps the layout", {
	expr <- cbind(c1=c(1, 2, 4, 8, 16, 100, 200, 300, 500),
				  c2=c(5, 5, 5, 5, 5, 1, 2, 3, 4))
	rownames(expr) <- paste0("g", 1:9)
	obj <- opsrun_obj(expr, chr=c(rep("chr1", 5), rep("chr2", 4)))
	res <- infercnv:::smooth_by_chromosome(obj, window_length=3)
	expect_equal(dim(res@expr.data), dim(expr))
	expect_equal(dimnames(res@expr.data), dimnames(obj@expr.data))
	expect_equal(res@expr.data[1:5, "c1"], setNames(c(4/3, 9/4, 4.5, 9, 40/3), paste0("g", 1:5)))
	## chr2 for cell c1: 100,200,300,500 with window 3
	## end1 = (2*100 + 200)/3 ; (100 + 2*200 + 300)/4 ; (200 + 600 + 500)/4 ; end4 = (300 + 2*500)/3
	expect_equal(unname(res@expr.data[6:9, "c1"]), c(400/3, 800/4, 1300/4, 1300/3))
	expect_equal(unname(res@expr.data[1:5, "c2"]), rep(5, 5))
	## everything else in the object is untouched
	expect_identical(res@gene_order, obj@gene_order)
	expect_identical(res@count.data, obj@count.data)
})

test_that("smooth_by_chromosome equals the reference per chromosome and cell, and ignores smooth_ends", {
	chr <- c(rep("chr1", 12), rep("chr2", 3), rep("chr3", 1), rep("chr4", 9))
	expr <- sapply(1:4, function(j) sin(seq_along(chr) * j) * 5 + j)
	obj <- opsrun_obj(expr, chr=chr)
	expected <- expr
	for (ch in unique(chr)) {
		i <- which(chr == ch)
		if (length(i) > 1) {
			for (j in 1:4) expected[i, j] <- opsrun_ref_pyramid(expr[i, j], 5)
		}
	}
	res <- infercnv:::smooth_by_chromosome(obj, window_length=5)
	expect_equal(unname(res@expr.data), unname(expected))
	## the single-gene chromosome is left as is
	expect_equal(res@expr.data[16, ], obj@expr.data[16, ])
	res2 <- infercnv:::smooth_by_chromosome(obj, window_length=5, smooth_ends=FALSE)
	expect_equal(res2@expr.data, res@expr.data)
})

test_that("smooth_by_chromosome with window_length < 2 leaves the data unchanged", {
	expr <- matrix(c(3, 1, 4, 1, 5, 9, 2, 6), ncol=2)
	obj <- opsrun_obj(expr)
	expect_equal(infercnv:::smooth_by_chromosome(obj, window_length=1)@expr.data, obj@expr.data)
})

test_that("smooth_by_chromosome mirrors the operation on the hspike object", {
	expr <- matrix(c(1, 2, 4, 8, 16, 3, 1, 4, 1, 5), ncol=2)
	hs <- opsrun_obj(matrix(c(10, 0, 10, 0, 10, 1, 2, 3, 4, 5), ncol=2))
	obj <- opsrun_with_hspike(opsrun_obj(expr), hs)
	res <- infercnv:::smooth_by_chromosome(obj, window_length=3)
	expect_equal(res@.hspike@expr.data, infercnv:::smooth_by_chromosome(hs, window_length=3)@expr.data)
	expect_equal(res@expr.data[, 1], setNames(c(4/3, 9/4, 4.5, 9, 40/3), paste0("g", 1:5)))
})

test_that("smooth_by_chromosome_runmeans is a running mean per chromosome (equal weights, truncated ends)", {
	chr <- c(rep("chr1", 7), rep("chr2", 4), "chr3")
	expr <- cbind(c1=c(1, 2, 4, 8, 16, 32, 64, 2, 4, 6, 8, 99),
				  c2=c(0, 0, 0, 9, 0, 0, 0, 1, 1, 1, 1, 7))
	obj <- opsrun_obj(expr, chr=chr)
	res <- infercnv:::smooth_by_chromosome_runmeans(obj, window_length=3)
	## hand values for chr1 / c1: 1.5, 7/3, 14/3, 28/3, 56/3, 112/3, 48
	expect_equal(unname(res@expr.data[1:7, "c1"]), c(1.5, 7/3, 14/3, 28/3, 56/3, 112/3, 48))
	expect_equal(unname(res@expr.data[8:11, "c1"]), opsrun_ref_runmean(c(2, 4, 6, 8), 3))
	## spike in c2 spreads over 3 genes and chromosomes do not leak into each other
	expect_equal(unname(res@expr.data[1:7, "c2"]), c(0, 0, 3, 3, 3, 0, 0))
	expect_equal(unname(res@expr.data[8:11, "c2"]), rep(1, 4))
	## the one-gene chromosome stays
	expect_equal(unname(res@expr.data[12, ]), c(99, 7))
	expect_equal(dimnames(res@expr.data), dimnames(obj@expr.data))
})

test_that("smooth_by_chromosome_runmeans matches the reference running mean for a longer window", {
	x <- sin(seq_len(30)) * 3 + seq_len(30) %% 4
	obj <- opsrun_obj(matrix(x, ncol=1))
	res <- infercnv:::smooth_by_chromosome_runmeans(obj, window_length=7)
	expect_equal(unname(res@expr.data[, 1]), opsrun_ref_runmean(x, 7))
})

test_that("smooth_by_chromosome_runmeans mirrors the hspike", {
	hs <- opsrun_obj(matrix(c(1, 2, 4, 8, 16), ncol=1))
	obj <- opsrun_with_hspike(opsrun_obj(matrix(c(5, 6, 7, 8, 9), ncol=1)), hs)
	res <- infercnv:::smooth_by_chromosome_runmeans(obj, window_length=3)
	expect_equal(unname(res@.hspike@expr.data[, 1]), opsrun_ref_runmean(c(1, 2, 4, 8, 16), 3))
})

test_that("smooth_by_chromosome_coordinates keeps constants and stays within the input range", {
	n <- 40
	start <- seq(1000, by=2000, length.out=n)
	x <- sin(seq_len(n)) * 3 + 10
	expr <- cbind(const=rep(2.5, n), var=x)
	obj <- opsrun_obj(expr, start=start)
	res <- infercnv:::smooth_by_chromosome_coordinates(obj, window_length=20000)
	expect_equal(dim(res@expr.data), dim(expr))
	expect_equal(dimnames(res@expr.data), dimnames(obj@expr.data))
	expect_equal(unname(res@expr.data[, "const"]), rep(2.5, n))
	expect_true(all(res@expr.data[, "var"] >= min(x) - 1e-12 & res@expr.data[, "var"] <= max(x) + 1e-12))
	## smoothing reduces the spread
	expect_lt(var(res@expr.data[, "var"]), var(x))
})

test_that("smooth_by_chromosome_coordinates smooths each chromosome separately", {
	n <- 12
	start <- rep(seq(1000, by=3000, length.out=6), 2)
	expr <- matrix(c(rep(1, 6), rep(100, 6)), ncol=1)
	obj <- opsrun_obj(expr, chr=rep(c("chr1", "chr2"), each=6), start=start)
	res <- infercnv:::smooth_by_chromosome_coordinates(obj, window_length=10000)
	expect_equal(unname(res@expr.data[1:6, 1]), rep(1, 6))
	expect_equal(unname(res@expr.data[7:12, 1]), rep(100, 6))
})

test_that("smooth_by_chromosome_coordinates leaves a one-gene chromosome alone", {
	expr <- matrix(c(1, 5, 9, 4), ncol=1)
	obj <- opsrun_obj(expr, chr=c("chr1", "chr1", "chr1", "chr2"), start=c(1000, 2000, 3000, 1000))
	res <- infercnv:::smooth_by_chromosome_coordinates(obj, window_length=10000)
	expect_equal(unname(res@expr.data[4, 1]), 4)
})

# pins the current output
## the first value is derived by hand: gene 1 has the window (-1250, 3750) so genes 1-3 with weights
## 1, 0.6, 0.2; the code then pads one more gene on the right with weight 0.1:
## (1*1 + 0.6*2 + 0.2*4 + 0.1*8) / 1.9 = 2
test_that("smooth_by_chromosome_coordinates output for a small example", {
	start <- c(1000, 2000, 3000, 4000, 5000)
	obj <- opsrun_obj(matrix(c(1, 2, 4, 8, 16), ncol=1), start=start)
	res <- infercnv:::smooth_by_chromosome_coordinates(obj, window_length=2500)
	expect_equal(unname(res@expr.data[, 1]),
				 c(2, 3.28, 5.15384615385, 8.2, 11.47368421053), tolerance=1e-9)
})

## ---------------------------------------------------------------------------------------------
## simple per-step transforms
## ---------------------------------------------------------------------------------------------

test_that("normalize_counts_by_seq_depth scales every cell to the median library size", {
	expr <- matrix(c(1, 3,   2, 6,   5, 5), nrow=2)
	## library sizes 4, 8, 10; median 8
	res <- infercnv:::normalize_counts_by_seq_depth(opsrun_obj(expr))
	expect_equal(unname(res@expr.data), matrix(c(2, 6,   2, 6,   4, 4), nrow=2))
	expect_equal(unname(colSums(res@expr.data)), c(8, 8, 8))
	## count.data is not changed
	expect_equal(unname(res@count.data), expr)
	## explicit factor
	res2 <- infercnv:::normalize_counts_by_seq_depth(opsrun_obj(expr), normalize_factor=100)
	expect_equal(unname(res2@expr.data), matrix(c(25, 75,   25, 75,   50, 50), nrow=2))
	expect_equal(dimnames(res@expr.data), dimnames(opsrun_obj(expr)@expr.data))
})

test_that(".normalize_data_matrix_by_seq_depth works on a plain matrix", {
	m <- matrix(c(2, 2, 0, 6, 1, 1), nrow=2, dimnames=list(c("a", "b"), c("x", "y", "z")))
	## sums 4, 6, 2; median 4
	expect_equal(infercnv:::.normalize_data_matrix_by_seq_depth(m),
				 matrix(c(2, 2, 0, 4, 2, 2), nrow=2, dimnames=dimnames(m)))
	expect_equal(infercnv:::.normalize_data_matrix_by_seq_depth(m, normalize_factor=1),
				 matrix(c(.5, .5, 0, 1, .5, .5), nrow=2, dimnames=dimnames(m)))
})

test_that("log2xplus1, invert_log2xplus1 and invert_log2 follow their definitions", {
	expr <- matrix(c(0, 1, 3, 7, 15, 255), nrow=2)
	obj <- opsrun_obj(expr)
	lg <- infercnv:::log2xplus1(obj)
	expect_equal(unname(lg@expr.data), matrix(c(0, 1, 2, 3, 4, 8), nrow=2))
	expect_equal(unname(infercnv:::invert_log2xplus1(lg)@expr.data), expr)
	expect_equal(unname(infercnv:::invert_log2(lg)@expr.data), matrix(c(1, 2, 4, 8, 16, 256), nrow=2))
	expect_equal(unname(lg@count.data), expr)
})

test_that("add_pseudocount and anscombe_transform follow their definitions", {
	obj <- opsrun_obj(matrix(c(0, 1, 5, 9.625), nrow=2))
	expect_equal(unname(infercnv:::add_pseudocount(obj, 0.5)@expr.data), matrix(c(0.5, 1.5, 5.5, 10.125), nrow=2))
	expect_equal(unname(infercnv:::add_pseudocount(obj, 0)@expr.data), unname(obj@expr.data))
	## 2*sqrt(x + 3/8): 0 -> 2*sqrt(.375) ; 1 -> 2*sqrt(1.375) ; 5 -> 2*sqrt(5.375) ; 9.625 -> 2*sqrt(10) 
	expect_equal(unname(infercnv:::anscombe_transform(obj)@expr.data),
				 matrix(c(2 * sqrt(0.375), 2 * sqrt(1.375), 2 * sqrt(5.375), 2 * sqrt(10)), nrow=2))
})

test_that("scale_infercnv_expr z-scores every gene across cells", {
	expr <- matrix(c(1, 10,  2, 20,  3, 30), nrow=2)
	res <- infercnv:::scale_infercnv_expr(opsrun_obj(expr))
	expect_equal(matrix(as.numeric(res@expr.data), nrow=2), matrix(c(-1, -1, 0, 0, 1, 1), nrow=2))
	expect_equal(dimnames(res@expr.data), dimnames(opsrun_obj(expr)@expr.data))
})

test_that("cross_cell_normalize matches upper quartiles across cells", {
	expr <- cbind(c(1, 2, 3, 4, 5), c(2, 4, 6, 8, 10))
	## upper quartiles 4 and 8, mean 6 -> multiply cells by 1.5 and 0.75
	res <- infercnv:::cross_cell_normalize(opsrun_obj(expr))
	expect_equal(unname(res@expr.data), cbind(c(1, 2, 3, 4, 5) * 1.5, c(1, 2, 3, 4, 5) * 1.5))
})

test_that("mean_center_gene_expr subtracts each gene's mean, ignoring NA", {
	expr <- matrix(c(1, 10,  2, 20,  6, NA), nrow=2)
	res <- infercnv:::mean_center_gene_expr(opsrun_obj(expr))
	expect_equal(unname(res@expr.data), matrix(c(-2, -5,  -1, 5,  3, NA), nrow=2))
})

test_that("apply_max_threshold_bounds clips values at +/- threshold", {
	expr <- matrix(c(-5, -1, 0, 1, 2.5, 3, 3.01, 10), nrow=2)
	res <- infercnv:::apply_max_threshold_bounds(opsrun_obj(expr), 3)
	expect_equal(unname(res@expr.data), matrix(c(-3, -1, 0, 1, 2.5, 3, 3, 3), nrow=2))
})

test_that("get_average_bounds averages the per-cell minimum and maximum", {
	expr <- cbind(c(1, 2, 3), c(10, 20, 30))
	expect_equal(infercnv:::.get_average_bounds(expr), c(5.5, 16.5))
	expect_equal(infercnv:::get_average_bounds(opsrun_obj(expr)), c(5.5, 16.5))
	## NA ignored
	expect_equal(infercnv:::.get_average_bounds(cbind(c(1, NA, 3), c(10, 20, 30))), c(5.5, 16.5))
})

test_that("transform_to_reference_based_Zscores centres on the reference and floors the sd at the mean", {
	## reference = cells 1,2 ; gene 1 ref (2,4): mean 3, sd 1.41 -> floored to 3 ; gene 2 ref (10,10): mean 10, sd 0 -> 10
	expr <- rbind(c(2, 4, 3, 9, 0),
				  c(10, 10, 11, 20, 0))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	res <- infercnv:::transform_to_reference_based_Zscores(obj)
	expect_equal(unname(res@expr.data), rbind(c(-1/3, 1/3, 0, 2, -1),
											  c(0, 0, 0.1, 1, -1)))
	## sd above the mean is used as is: ref (0, 10): mean 5, sd 7.07 -> (x - 5)/sd
	obj2 <- opsrun_obj(rbind(c(0, 10, 5, 15), c(1, 1, 1, 1)), ref=list(r=c(1, 2)))
	res2 <- infercnv:::transform_to_reference_based_Zscores(obj2)
	expect_equal(unname(res2@expr.data[1, ]), (c(0, 10, 5, 15) - 5) / sd(c(0, 10)))
	## second gene: ref mean 1, sd 0 -> floored to 1
	expect_equal(unname(res2@expr.data[2, ]), rep(0, 4))
})

test_that("get_reference_grouped_cell_indices concatenates all reference groups", {
	obj <- opsrun_obj(matrix(1, 2, 6), ref=list(a=c(1, 2), b=c(5)))
	expect_equal(unname(infercnv:::get_reference_grouped_cell_indices(obj)), c(1, 2, 5))
	expect_equal(names(infercnv:::get_reference_grouped_cell_indices(obj)), c("a1", "a2", "b"))
	expect_null(infercnv:::get_reference_grouped_cell_indices(opsrun_obj(matrix(1, 2, 3))))
})

test_that("simple transforms are mirrored on the hspike object", {
	ex <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow=3)
	hx <- matrix(c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11), nrow=3)
	obj <- opsrun_with_hspike(opsrun_obj(ex, ref=list(r=c(1, 2))), opsrun_obj(hx, ref=list(r=c(1, 2))))
	hs <- opsrun_obj(hx, ref=list(r=c(1, 2)))
	ops <- list(infercnv:::log2xplus1,
				infercnv:::invert_log2xplus1,
				infercnv:::invert_log2,
				function(o) infercnv:::add_pseudocount(o, 2),
				infercnv:::anscombe_transform,
				infercnv:::scale_infercnv_expr,
				infercnv:::cross_cell_normalize,
				infercnv:::mean_center_gene_expr,
				function(o) infercnv:::apply_max_threshold_bounds(o, 5),
				infercnv:::transform_to_reference_based_Zscores,
				function(o) infercnv:::center_cell_expr_across_chromosome(o, "median"),
				function(o) infercnv:::subtract_ref_expr_from_obs(o),
				function(o) infercnv:::remove_outliers_norm(o))
	for (i in seq_along(ops)) {
		expect_equal(ops[[i]](obj)@.hspike@expr.data, ops[[i]](hs)@expr.data, info=paste("op", i))
	}
})

## ---------------------------------------------------------------------------------------------
## gene filters
## ---------------------------------------------------------------------------------------------

test_that("require_above_min_mean_expr_cutoff drops genes whose mean is strictly below the cutoff", {
	expr <- rbind(g1=c(0, 0, 0), g2=c(1, 1, 1), g3=c(2, 0, 1), g4=c(5, 5, 5))
	obj <- opsrun_obj(expr, chr=c("chr1", "chr1", "chr2", "chr2"))
	res <- infercnv:::require_above_min_mean_expr_cutoff(obj, 1)
	expect_equal(rownames(res@expr.data), c("g2", "g3", "g4"))
	expect_equal(rownames(res@count.data), c("g2", "g3", "g4"))
	expect_equal(rownames(res@gene_order), c("g2", "g3", "g4"))
	expect_equal(levels(res@gene_order$chr), c("chr1", "chr2"))
	res2 <- infercnv:::require_above_min_mean_expr_cutoff(obj, 1.5)
	expect_equal(rownames(res2@expr.data), "g4")
	## dropped chromosomes lose their factor level
	expect_equal(levels(res2@gene_order$chr), "chr2")
	## nothing to remove -> the object is returned as is
	expect_identical(infercnv:::require_above_min_mean_expr_cutoff(obj, 0), obj)
	expect_equal(unname(infercnv:::.below_min_mean_expr_cutoff(expr, 1.5)), 1:3)
})

test_that("require_above_min_cells_ref keeps genes with enough cells above zero (NA not counted)", {
	expr <- rbind(g1=c(0, 0, 1, 0), g2=c(1, 1, 0, 0), g3=c(1, NA, 2, 3), g4=c(NA, NA, 5, 0), g5=c(1, 1, 1, 1))
	obj <- opsrun_obj(expr, chr=c("chr1", "chr1", "chr1", "chr2", "chr2"))
	res <- infercnv:::require_above_min_cells_ref(obj, min_cells_per_gene=2)
	expect_equal(rownames(res@expr.data), c("g2", "g3", "g5"))
	expect_equal(rownames(res@gene_order), c("g2", "g3", "g5"))
	expect_equal(levels(res@gene_order$chr), c("chr1", "chr2"))
	res3 <- infercnv:::require_above_min_cells_ref(obj, min_cells_per_gene=3)
	expect_equal(rownames(res3@expr.data), c("g3", "g5"))
	## exactly the minimum number of cells is enough
	expect_equal(rownames(infercnv:::require_above_min_cells_ref(obj, 1)@expr.data), paste0("g", 1:5))
	expect_identical(infercnv:::require_above_min_cells_ref(obj, 0), obj)
})

test_that("require_above_min_cells_ref stops when every gene would be removed", {
	obj <- opsrun_obj(matrix(c(1, 0, 0, 1, 0, 0), nrow=3))
	expect_error(infercnv:::require_above_min_cells_ref(obj, 5), "998")
})

## ---------------------------------------------------------------------------------------------
## subtracting the reference
## ---------------------------------------------------------------------------------------------

test_that(".get_normal_gene_mean_bounds gives per gene and group means (plain and inverse-log)", {
	expr <- rbind(g1=c(1, 3, 10), g2=c(0, 0, 4))
	groups <- list(r1=c(1, 2), r2=3)
	res <- infercnv:::.get_normal_gene_mean_bounds(expr, groups)
	expect_equal(rownames(res), c("g1", "g2"))
	expect_equal(colnames(res), c("r1", "r2"))
	expect_equal(res[["r1"]], c(2, 0))
	expect_equal(res[["r2"]], c(10, 4))
	## inv_log: log2(mean(2^x - 1) + 1) ; g1/r1: 2^1-1 = 1, 2^3-1 = 7, mean 4 -> log2(5)
	res_il <- infercnv:::.get_normal_gene_mean_bounds(expr, groups, inv_log=TRUE)
	expect_equal(res_il[["r1"]], c(log2(5), 0))
	expect_equal(res_il[["r2"]], c(10, 4))
})

test_that(".subtract_expr subtracts group means with and without bounds", {
	expr <- matrix(c(10, 20,   4, 4,   0, 30,   6, 12), nrow=2,
				   dimnames=list(c("g1", "g2"), paste0("c", 1:4)))
	means <- data.frame(r1=c(5, 8), r2=c(9, 8), row.names=c("g1", "g2"))
	## without bounds: minus the mean of the group means (7, 8)
	expect_equal(infercnv:::.subtract_expr(expr, means, use_bounds=FALSE),
				 matrix(c(3, 12,  -3, -4,  -7, 22,  -1, 4), nrow=2, dimnames=dimnames(expr)))
	## with bounds: only the part outside [min, max] of the group means remains; g1 range 5..9, g2 range 8..8
	expect_equal(infercnv:::.subtract_expr(expr, means, use_bounds=TRUE),
				 matrix(c(1, 12,  -1, -4,  -5, 22,  0, 4), nrow=2, dimnames=dimnames(expr)))
})

test_that("subtract_ref_expr_from_obs handles inv_log, use_bounds=FALSE and a missing reference", {
	## cells 1,2 are reference
	expr <- rbind(g1=c(1, 3, 6, 0), g2=c(2, 2, 2, 5))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	res <- infercnv:::subtract_ref_expr_from_obs(obj, use_bounds=FALSE)
	expect_equal(unname(res@expr.data), rbind(c(-1, 1, 4, -2), c(0, 0, 0, 3)))
	expect_equal(rownames(res@expr.data), c("g1", "g2"))
	## inv_log: ref mean on the linear scale log2(mean(2^c(1,3) - 1) + 1) = log2(5)
	res_il <- infercnv:::subtract_ref_expr_from_obs(obj, inv_log=TRUE, use_bounds=FALSE)
	expect_equal(unname(res_il@expr.data[1, ]), c(1, 3, 6, 0) - log2(5))
	expect_equal(unname(res_il@expr.data[2, ]), c(0, 0, 0, 3))
	## no reference cells: mean of all cells is the proxy
	obj_noref <- opsrun_obj(expr, obs=list(a=1:2, b=3:4))
	res_nr <- infercnv:::subtract_ref_expr_from_obs(obj_noref, use_bounds=FALSE)
	expect_equal(unname(res_nr@expr.data), rbind(c(1, 3, 6, 0) - mean(c(1, 3, 6, 0)), c(2, 2, 2, 5) - 2.75))
	## with bounds and no reference, a one-valued range [mean, mean] leaves cells above or below it
	res_nrb <- infercnv:::subtract_ref_expr_from_obs(obj_noref, use_bounds=TRUE)
	expect_equal(unname(res_nrb@expr.data), unname(res_nr@expr.data))
})

test_that("subtract_ref_expr_from_obs with several reference groups uses the min and max of group means", {
	expr <- rbind(g1=c(0, 10, 4, 5, 20, -3))
	obj <- opsrun_obj(expr, ref=list(a=c(1, 2), b=c(3)))
	## group means 5 and 4 -> band [4, 5]; ref-free cells 4,5,6 are 5, 20, -3
	res <- infercnv:::subtract_ref_expr_from_obs(obj, use_bounds=TRUE)
	expect_equal(unname(res@expr.data[1, ]), c(-4, 5, 0, 0, 15, -7))
})

test_that("split_references splits reference cells into the requested number of clusters", {
	expr <- rbind(g1=c(0, 100, 1, 101, 5),
				  g2=c(0, 100, 0, 100, 5),
				  g3=c(1, 99, 0, 100, 5))
	colnames(expr) <- paste0("c", 1:5)
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2, 3, 4)), obs=list(t=5))
	res <- infercnv:::split_references(obj, num_groups=2, hclust_method="complete")
	expect_equal(names(res@reference_grouped_cell_indices), c("refgrp-1", "refgrp-2"))
	expect_equal(res@reference_grouped_cell_indices[["refgrp-1"]], c(1L, 3L))
	expect_equal(res@reference_grouped_cell_indices[["refgrp-2"]], c(2L, 4L))
	## expression and observations are unchanged
	expect_identical(res@expr.data, obj@expr.data)
	expect_identical(res@observation_grouped_cell_indices, obj@observation_grouped_cell_indices)
	## one group keeps all reference cells together
	res1 <- infercnv:::split_references(obj, num_groups=1, hclust_method="complete")
	expect_equal(names(res1@reference_grouped_cell_indices), "refgrp-1")
	expect_equal(res1@reference_grouped_cell_indices[["refgrp-1"]], 1:4)
})

## ---------------------------------------------------------------------------------------------
## centering, chromosome ends, outliers
## ---------------------------------------------------------------------------------------------

test_that(".center_columns centres each cell on its median or mean and ignores NA", {
	m <- cbind(a=c(1, 2, 10), b=c(4, NA, 8), c=c(5, 5, 5))
	rownames(m) <- c("g1", "g2", "g3")
	med <- infercnv:::.center_columns(m, "median")
	expect_equal(unname(med), cbind(c(-1, 0, 8), c(-2, NA, 2), c(0, 0, 0)))
	expect_equal(dimnames(med), dimnames(m))
	mn <- infercnv:::.center_columns(m, "mean")
	expect_equal(unname(mn), cbind(c(1 - 13/3, 2 - 13/3, 10 - 13/3), c(-2, NA, 2), c(0, 0, 0)))
})

test_that("center_cell_expr_across_chromosome centres whole cells, not chromosomes", {
	expr <- cbind(c(1, 2, 3, 10, 20, 30), c(0, 0, 0, 5, 5, 5))
	obj <- opsrun_obj(expr, chr=rep(c("chr1", "chr2"), each=3))
	res <- infercnv:::center_cell_expr_across_chromosome(obj, method="median")
	## medians 6.5 and 2.5
	expect_equal(unname(res@expr.data), cbind(expr[, 1] - 6.5, expr[, 2] - 2.5))
	res_mean <- infercnv:::center_cell_expr_across_chromosome(obj)
	expect_equal(unname(res_mean@expr.data), cbind(expr[, 1] - 11, expr[, 2] - 2.5))
})

test_that("remove_genes_at_ends_of_chromosomes drops (window_length - 1)/2 genes at both chromosome ends", {
	chr <- c(rep("chr1", 10), rep("chr2", 5), rep("chr3", 2))
	obj <- opsrun_obj(matrix(seq_along(chr), ncol=1), chr=chr)
	res <- infercnv:::remove_genes_at_ends_of_chromosomes(obj, window_length=7)
	## chr1 loses 3 + 3 genes, chr2 (5 < 2*3) loses floor(5/3) = 1 + 1 genes, chr3 (< 3 genes) none
	expect_equal(rownames(res@expr.data), paste0("g", c(4:7, 12:14, 16:17)))
	expect_equal(rownames(res@gene_order), rownames(res@expr.data))
	expect_equal(rownames(res@count.data), rownames(res@expr.data))
	expect_equal(levels(res@gene_order$chr), c("chr1", "chr2", "chr3"))
})

test_that("remove_genes_at_ends_of_chromosomes stops when no gene would be removed", {
	obj <- opsrun_obj(matrix(1:4, ncol=1), chr=c("chr1", "chr1", "chr2", "chr2"))
	expect_error(infercnv:::remove_genes_at_ends_of_chromosomes(obj, window_length=7), "1234")
})

test_that("remove_genes_at_ends_of_chromosomes mirrors the removal on the hspike", {
	chr <- rep("chr1", 10)
	obj <- opsrun_with_hspike(opsrun_obj(matrix(1:10, ncol=1), chr=chr),
							  opsrun_obj(matrix(11:20, ncol=1), chr=chr))
	res <- infercnv:::remove_genes_at_ends_of_chromosomes(obj, window_length=5 + 2)
	expect_equal(nrow(res@.hspike@expr.data), 4)
	expect_equal(unname(res@.hspike@expr.data[, 1]), 14:17)
})

test_that("remove_outliers_norm with average_bound clips at the mean per-cell min and max", {
	## per cell: min 1/max 5 and min 0/max 40 -> bounds 0.5 and 22.5
	expr <- cbind(c(1, 2, 3, 4, 5), c(0, 10, 20, 30, 40))
	res <- infercnv:::.remove_outliers_norm(expr)
	expect_equal(res, cbind(c(1, 2, 3, 4, 5), c(0.5, 10, 20, 22.5, 22.5)))
	obj <- opsrun_obj(expr)
	expect_equal(unname(infercnv:::remove_outliers_norm(obj)@expr.data), unname(res))
})

test_that(".remove_outliers_norm uses fixed bounds only when both are given", {
	expr <- cbind(c(1, 2, 3, 4, 5), c(0, 10, 20, 30, 40))
	expect_equal(infercnv:::.remove_outliers_norm(expr, lower_bound=2, upper_bound=25),
				 cbind(c(2, 2, 3, 4, 5), c(2, 10, 20, 25, 25)))
	## a single bound is ignored, the method decides
	expect_equal(infercnv:::.remove_outliers_norm(expr, lower_bound=2),
				 infercnv:::.remove_outliers_norm(expr))
	expect_equal(infercnv:::.remove_outliers_norm(expr, out_method=NA, lower_bound=2, upper_bound=25),
				 cbind(c(2, 2, 3, 4, 5), c(2, 10, 20, 25, 25)))
})

test_that(".remove_outliers_norm stops on bad input", {
	expr <- cbind(c(1, 2, 3), c(4, 5, 6))
	expect_error(infercnv:::.remove_outliers_norm(expr, out_method="something_else"), "991")
	expect_error(infercnv:::.remove_outliers_norm(expr, out_method=NA), "992")
	expect_error(infercnv:::.remove_outliers_norm(NULL), "something is wrong with the data")
	expect_error(infercnv:::.remove_outliers_norm(matrix(numeric(0), 0, 2)), "something is wrong with the data")
})

## ---------------------------------------------------------------------------------------------
## noise removal
## ---------------------------------------------------------------------------------------------

test_that(".clear_noise sets values strictly inside (centre - threshold, centre + threshold) to the centre", {
	x <- matrix(c(0.5, 0.9, 1, 1.9, 2, 2.1, 3, 3.01), nrow=2)
	## centre 2, threshold 1 -> (1, 3) exclusive
	expect_equal(infercnv:::.clear_noise(x, threshold=1, center_pos=2),
				 matrix(c(0.5, 0.9, 1, 2, 2, 2, 3, 3.01), nrow=2))
})

test_that("clear_noise centres on the reference mean and leaves values at the boundary alone", {
	## reference cells 1,2 : mean 3 ; threshold 1 -> (2, 4) exclusive
	expr <- cbind(c(1, 2, 3), c(3, 4, 5), c(2, 2.5, 4), c(3.9, 6, 0))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	res <- infercnv:::clear_noise(obj, threshold=1)
	expect_equal(unname(res@expr.data), cbind(c(1, 2, 3), c(3, 4, 5), c(2, 3, 4), c(3, 6, 0)))
	## threshold 0 returns the object untouched
	expect_identical(infercnv:::clear_noise(obj, threshold=0), obj)
})

test_that("clear_noise without reference cells uses the mean of all data", {
	expr <- cbind(c(0, 10), c(5, 5))
	obj <- opsrun_obj(expr, obs=list(a=1:2))
	## mean 5, threshold 2 -> values in (3, 7): the two 5s are already at the centre
	res <- infercnv:::clear_noise(obj, threshold=2)
	expect_equal(unname(res@expr.data), expr)
	res2 <- infercnv:::clear_noise(opsrun_obj(cbind(c(0, 4), c(6, 10)), obs=list(a=1:2)), threshold=2)
	## mean 5: 4 and 6 are within (3, 7)
	expect_equal(unname(res2@expr.data), cbind(c(0, 5), c(5, 10)))
})

test_that("clear_noise with noise_logistic depresses values close to the centre smoothly", {
	## reference mean is 0
	expr <- cbind(c(-1, 1, 0), c(1, -1, 0), c(0.5, -0.5, 1.5))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	res <- infercnv:::clear_noise(obj, threshold=0.5, noise_logistic=TRUE)
	p <- function(v) 1 / (1 + exp(-20 * (abs(v) - 0.5)))
	## new = centre +/- p(|x - centre|) * |x - centre| with a logistic of slope 20 and midpoint = threshold
	exp_vals <- sign(expr) * p(expr) * abs(expr)
	expect_equal(unname(res@expr.data), exp_vals)
	## hand-worked: exactly at the midpoint the logistic is 1/2
	expect_equal(unname(res@expr.data[1, 3]), 0.25)
	expect_equal(unname(res@expr.data[2, 3]), -0.25)
	## far above the midpoint the value is almost unchanged, the centre stays put
	expect_equal(unname(res@expr.data[3, 3]), 1.5 * p(1.5))
	expect_equal(unname(res@expr.data[3, 1]), 0)
})

test_that("clear_noise_via_ref_mean_sd removes values within mean +/- sd_amplifier * mean per-cell sd", {
	## ref cells: (1,2,3) sd 1 and (2,4,6) sd 2 -> mean sd 1.5 ; ref mean 3
	expr <- cbind(c(1, 2, 3), c(2, 4, 6), c(0.8, 5.2, 5.3))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	## amplifier 1: band (1.5, 4.5) exclusive
	res <- infercnv:::clear_noise_via_ref_mean_sd(obj, sd_amplifier=1)
	expect_equal(unname(res@expr.data), cbind(c(1, 3, 3), c(3, 3, 6), c(0.8, 5.2, 5.3)))
	## amplifier 2: band (0, 6) exclusive
	res2 <- infercnv:::clear_noise_via_ref_mean_sd(obj, sd_amplifier=2)
	expect_equal(unname(res2@expr.data), cbind(c(3, 3, 3), c(3, 3, 6), c(3, 3, 3)))
	## default amplifier is 1.5: band 3 +/- 2.25 = (0.75, 5.25)
	res3 <- infercnv:::clear_noise_via_ref_mean_sd(obj)
	expect_equal(unname(res3@expr.data), cbind(c(3, 3, 3), c(3, 3, 6), c(3, 3, 5.3)))
})

test_that("clear_noise_via_ref_mean_sd uses all observations when there is no reference", {
	expr <- cbind(c(0, 2, 4), c(2, 4, 6))
	obj <- opsrun_obj(expr, obs=list(a=1:2))
	## mean 3, per-cell sd 2 each -> band 3 +/- 2 -> (1, 5): 2, 4, 2, 4 -> 3
	res <- infercnv:::clear_noise_via_ref_mean_sd(obj, sd_amplifier=1)
	expect_equal(unname(res@expr.data), cbind(c(0, 3, 3), c(3, 3, 6)))
})

test_that("clear_noise_via_ref_mean_sd with noise_logistic uses the sd band as the midpoint", {
	expr <- cbind(c(1, 2, 3), c(2, 4, 6), c(0, 3, 6))
	obj <- opsrun_obj(expr, ref=list(r=c(1, 2)))
	res <- infercnv:::clear_noise_via_ref_mean_sd(obj, sd_amplifier=1, noise_logistic=TRUE)
	centre <- 3
	delta <- 1.5
	expected <- centre + sign(expr - centre) / (1 + exp(-20 * (abs(expr - centre) - delta))) * abs(expr - centre)
	expect_equal(unname(res@expr.data), expected)
})

test_that(".apply_logistic_val_adj with a given midpoint and slope", {
	expr <- matrix(c(2, 4, 6, 8), nrow=2)
	## centre 5, delta 1 slope 1: |x - 5| = 3, 1, 1, 3
	res <- infercnv:::.apply_logistic_val_adj(expr, expr_mean=5, delta_midpt=1, slope=1)
	l3 <- 1 / (1 + exp(-2))
	expect_equal(unname(res), matrix(c(5 - l3 * 3, 5 - 0.5 * 1, 5 + 0.5 * 1, 5 + l3 * 3), nrow=2))
	expect_equal(infercnv:::.logistic(1, 1, 5), 0.5)
	expect_equal(infercnv:::.logistic(0, 1, 1), 1 / (1 + exp(1)))
})

## ---------------------------------------------------------------------------------------------
## re-clustering after the subcluster definition
## ---------------------------------------------------------------------------------------------

test_that(".redo_hierarchical_clustering recomputes a tree for the cells of every group", {
	expr <- rbind(g1=c(0, 50, 1, 51, 5, 100),
				  g2=c(1, 52, 0, 50, 90, 110),
				  g3=c(2, 49, 1, 49, 4, 3))
	colnames(expr) <- paste0("c", 1:6)
	obj <- opsrun_obj(expr, obs=list(a=1:4, b=5:6))
	obj@tumor_subclusters <- list(subclusters=list(a=list(a_s1=c(c1=1L, c3=3L), a_s2=c(c4=4L, c2=2L)),
												   b=list(b_s1=c(c5=5L, c6=6L))))
	res <- infercnv:::.redo_hierarchical_clustering(obj, hclust_method="complete")
	expect_equal(sort(names(res@tumor_subclusters$hc)), c("a", "b"))
	ref_a <- stats::hclust(stats::dist(t(expr[, 1:4])), method="complete")
	expect_s3_class(res@tumor_subclusters$hc$a, "hclust")
	expect_equal(res@tumor_subclusters$hc$a$height, ref_a$height)
	expect_equal(res@tumor_subclusters$hc$a$merge, ref_a$merge)
	expect_equal(res@tumor_subclusters$hc$a$labels, paste0("c", 1:4))
	ref_b <- stats::hclust(stats::dist(t(expr[, 5:6])), method="complete")
	expect_equal(res@tumor_subclusters$hc$b$height, ref_b$height)
	## the subclusters themselves are untouched
	expect_identical(res@tumor_subclusters$subclusters, obj@tumor_subclusters$subclusters)
})

## ---------------------------------------------------------------------------------------------
## resume-mode bookkeeping used by run()
## ---------------------------------------------------------------------------------------------

test_that("compareNA is TRUE for equal vectors including NA positions", {
	expect_true(infercnv:::compareNA(c(1, NA, 3), c(1, NA, 3)))
	expect_true(infercnv:::compareNA("a", "a"))
	expect_true(infercnv:::compareNA(NULL, NULL))
	expect_true(infercnv:::compareNA(NA, NA))
	expect_false(infercnv:::compareNA(c(1, 2), c(1, 3)))
	expect_false(infercnv:::compareNA(c(1, 2), c(1, 2, 3)))
	expect_false(infercnv:::compareNA(NULL, 1))
	expect_false(infercnv:::compareNA(1, NULL))
})

test_that(".compare_args compares only the relevant arguments", {
	cur <- list(a=1, b="x", c=NA)
	expect_true(infercnv:::.compare_args(cur, c("a", "b"), list(a=1, b="x", c=5)))
	expect_false(infercnv:::.compare_args(cur, c("a", "b"), list(a=2, b="x", c=5)))
	expect_false(infercnv:::.compare_args(cur, c("a", "b"), list(a=1, b="y")))
	## an irrelevant difference is ignored
	expect_true(infercnv:::.compare_args(cur, "a", list(a=1, b="other", c=7)))
	## NA equals NA
	expect_true(infercnv:::.compare_args(cur, "c", list(a=9, c=NA)))
	## nothing relevant -> equal
	expect_true(infercnv:::.compare_args(cur, character(0), list(a=100)))
})

test_that(".compare_args is FALSE when a relevant argument is missing on one side", {
	expect_false(infercnv:::.compare_args(list(a=1, b=2), c("a", "b"), list(a=1)))
	expect_false(infercnv:::.compare_args(list(a=1), c("a", "b"), list(a=1, b=2)))
	## a relevant argument missing on both sides is not compared
	expect_true(infercnv:::.compare_args(list(a=1), c("a", "zzz"), list(a=1)))
})

opsrun_default_run_args <- function(...) {
	args <- list(out_dir="od", HMM=FALSE, HMM_type="i6", num_ref_groups=NULL, analysis_mode="subclusters",
				 tumor_subcluster_partition_method="leiden", smooth_method="pyramidinal",
				 max_centered_threshold=3, remove_genes_at_chr_ends=FALSE, prune_outliers=FALSE,
				 BayesMaxPNormal=0.5, mask_nonDE_genes=FALSE, denoise=FALSE, noise_filter=NA,
				 sd_amplifier=1.5, noise_logistic=FALSE)
	over <- list(...)
	for (n in names(over)) args[n] <- list(over[[n]])
	args
}

test_that(".get_relevant_args_list names the step files for a default run", {
	res <- do.call(infercnv:::.get_relevant_args_list, opsrun_default_run_args())
	expected <- c("01_incoming_data.infercnv_obj",
				  "02_reduced_by_cutoff.infercnv_obj",
				  "03_normalized_by_depth.infercnv_obj",
				  "04_logtransformed.infercnv_obj",
				  "05_scaled.infercnv_obj",
				  "",
				  "",
				  "08_remove_ref_avg_from_obs_logFC.infercnv_obj",
				  "09_apply_max_centered_expr_threshold.infercnv_obj",
				  "10_smoothed_by_chr.infercnv_obj",
				  "11_recentered_cells_by_chr.infercnv_obj",
				  "12_remove_ref_avg_from_obs_adjust.infercnv_obj",
				  "",
				  "14_invert_log_transform.infercnv_obj",
				  "15_tumor_subclusters.leiden.infercnv_obj",
				  "", "", "", "", "", "")
	expect_equal(res$expected_file_names, ifelse(expected == "", "", file.path("od", expected)))
	expect_equal(res$relevant_args[[1]], c("chr_exclude", "max_cells_per_group", "min_max_counts_per_cell", "counts_md5"))
	expect_equal(res$relevant_args[[2]], c("cutoff", "min_cells_per_gene"))
	expect_equal(res$relevant_args[[3]], "HMM")
	expect_equal(res$relevant_args[[10]], c("smooth_method", "window_length", "smooth_ends"))
	expect_equal(res$relevant_args[[15]][1:4], c("analysis_mode", "tumor_subcluster_partition_method", "z_score_filter", "k_nn"))
	expect_true(all(c("hclust_method", "cluster_by_groups", "leiden_resolution") %in% res$relevant_args[[15]]))
})

test_that(".get_relevant_args_list adds HMM, reference-split, subcluster-less, outlier, mask and denoise steps", {
	res <- do.call(infercnv:::.get_relevant_args_list,
				   opsrun_default_run_args(HMM=TRUE, analysis_mode="samples", num_ref_groups=2, denoise=TRUE,
										   mask_nonDE_genes=TRUE, prune_outliers=TRUE, remove_genes_at_chr_ends=TRUE))
	expected <- c("01_incoming_data.infercnv_obj",
				  "02_reduced_by_cutoff.infercnv_obj",
				  "03_normalized_by_depthHMMi6.infercnv_obj",
				  "04_logtransformedHMMi6.infercnv_obj",
				  "05_scaledHMMi6.infercnv_obj",
				  "06_split_HMMi6f_refs2.infercnv_obj",
				  "",
				  "08_remove_ref_avg_from_obs_logFCHMMi6.infercnv_obj",
				  "09_apply_max_centered_expr_thresholdHMMi6.infercnv_obj",
				  "10_smoothed_by_chrHMMi6.infercnv_obj",
				  "11_recentered_cells_by_chrHMMi6.infercnv_obj",
				  "12_remove_ref_avg_from_obs_adjustHMMi6.infercnv_obj",
				  "13_remove_gene_at_chr_endsHMMi6.infercnv_obj",
				  "14_invert_log_transformHMMi6.infercnv_obj",
				  "15_no_subclusteringHMMi6.infercnv_obj",
				  "16_remove_outlierHMMi6.infercnv_obj",
				  "17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj",
				  "18_HMM_pred.Bayes_NetHMMi6.hmm_mode-samples.mcmc_obj",
				  "19_HMM_pred.Bayes_NetHMMi6.hmm_mode-samples.Pnorm_0.5.infercnv_obj",
				  "20_HMM_pred.repr_intensitiesHMMi6.hmm_mode-samples.Pnorm_0.5.infercnv_obj",
				  "21_mask_nonDEHMMi6.infercnv_obj",
				  "22_denoiseHMMi6.NF_NA.SD_1.5.NL_FALSE.infercnv_obj")
	expect_equal(res$expected_file_names, ifelse(expected == "", "", file.path("od", expected)))
	expect_equal(res$relevant_args[[3]], c("HMM", "HMM_type", "sim_method", "hspike_aggregate_normals", "sim_foreground"))
	expect_equal(res$relevant_args[[6]], c("num_ref_groups", "hclust_method"))
	expect_equal(res$relevant_args[[13]], c("remove_genes_at_chr_ends", "window_length"))
	expect_equal(res$relevant_args[[16]], c("prune_outliers", "outlier_method_bound", "outlier_lower_bound", "outlier_upper_bound"))
	expect_equal(res$relevant_args[[21]], c("mask_nonDE_genes", "require_DE_all_normals", "test.use", "mask_nonDE_pval"))
	expect_equal(res$relevant_args[[22]], c("denoise", "noise_filter", "sd_amplifier", "noise_logistic"))
})

test_that(".get_relevant_args_list names the denoise file after the settings and omits max-threshold file when NA", {
	res <- do.call(infercnv:::.get_relevant_args_list,
				   opsrun_default_run_args(denoise=TRUE, noise_filter=0.2, sd_amplifier=2, noise_logistic=TRUE,
										   max_centered_threshold=NA, tumor_subcluster_partition_method="qnorm"))
	expect_equal(res$expected_file_names[[22]], file.path("od", "22_denoise.qnorm.NF_0.2.SD_2.NL_TRUE.infercnv_obj"))
	expect_equal(res$expected_file_names[[9]], "")
	expect_equal(res$expected_file_names[[15]], file.path("od", "15_tumor_subclusters.qnorm.infercnv_obj"))
})

test_that(".get_relevant_args_list random_trees moves the subcluster file to step 7", {
	res <- do.call(infercnv:::.get_relevant_args_list,
				   opsrun_default_run_args(tumor_subcluster_partition_method="random_trees"))
	expect_equal(res$expected_file_names[[7]], file.path("od", "07_tumor_subclusters.rand_trees.random_trees.infercnv_obj"))
	expect_equal(res$expected_file_names[[8]], file.path("od", "08_remove_ref_avg_from_obs_logFC.rand_trees.infercnv_obj"))
	expect_equal(res$expected_file_names[[15]], "")
	expect_true(all(c("hclust_method", "tumor_subcluster_pval", "cluster_by_groups") %in% res$relevant_args[[7]]))
})

## ---------------------------------------------------------------------------------------------
## median filtering (R/noise_reduction.R)
## ---------------------------------------------------------------------------------------------

## The windows below are chosen so that every window covers the whole block (blocks have at most 5 genes
## and 5 cells), so each entry becomes the median of its block, whatever the exact window reach is.

test_that(".median_filter replaces every entry by the block median when the window covers the block", {
	d <- matrix(c(1, 5, 3, 9, 2, 8), nrow=3)
	res <- infercnv:::.median_filter(d, window_size=7, half_window=3)
	expect_equal(res, matrix(4, 3, 2))
	## a constant matrix stays unchanged for any window
	k <- matrix(2.5, 12, 6)
	expect_equal(infercnv:::.median_filter(k, window_size=3, half_window=1), k)
	expect_equal(infercnv:::.median_filter(k, window_size=7, half_window=3), k)
})

test_that(".median_filter keeps dimensions and names", {
	d <- matrix(1:6, nrow=3, dimnames=list(c("a", "b", "c"), c("x", "y")))
	res <- infercnv:::.median_filter(d, 7, 3)
	expect_equal(dimnames(res), dimnames(d))
	expect_equal(unname(res), matrix(3.5, 3, 2))
})

test_that("apply_median_filtering filters tumor subclusters and reference groups separately, per chromosome", {
	## 6 genes (chr1: g1-g3, chr2: g4-g6) x 8 cells: refs = cells 1-2, tumor subclusters s1 = cells 3-5, s2 = cells 6-8
	expr <- matrix(0, 6, 8)
	expr[1:3, 1:2] <- c(1, 2, 3, 10, 20, 30)
	expr[4:6, 1:2] <- c(100, 200, 300, 4, 5, 6)
	expr[1:3, 3:5] <- c(1, 1, 1, 2, 2, 2, 3, 3, 9)
	expr[4:6, 3:5] <- 7
	expr[1:3, 6:8] <- 1:9
	expr[4:6, 6:8] <- 50
	obj <- opsrun_obj(expr, chr=rep(c("chr1", "chr2"), each=3), ref=list(r=1:2), obs=list(t=3:8))
	obj@tumor_subclusters <- list(subclusters=list(t=list(t_s1=3:5, t_s2=6:8)))
	res <- infercnv:::apply_median_filtering(obj, window_size=7)
	expected <- expr
	expected[1:3, 1:2] <- median(c(1, 2, 3, 10, 20, 30))
	expected[4:6, 1:2] <- median(c(100, 200, 300, 4, 5, 6))
	expected[1:3, 3:5] <- median(c(1, 1, 1, 2, 2, 2, 3, 3, 9))
	expected[4:6, 3:5] <- 7
	expected[1:3, 6:8] <- median(1:9)
	expected[4:6, 6:8] <- 50
	expect_equal(unname(res@expr.data), expected)
	expect_equal(dimnames(res@expr.data), dimnames(obj@expr.data))
	## nothing else changes
	expect_identical(res@tumor_subclusters, obj@tumor_subclusters)
	expect_identical(res@count.data, obj@count.data)
})

test_that("apply_median_filtering can restrict itself to observations or references", {
	expr <- rbind(c(1, 2, 3, 10, 20, 30),
				  c(4, 5, 6, 11, 21, 31),
				  c(7, 8, 9, 12, 22, 32))
	obj <- opsrun_obj(expr, ref=list(r=1:3), obs=list(t=4:6))
	obj@tumor_subclusters <- list(subclusters=list(t=list(t_s1=4:6)))
	only_obs <- infercnv:::apply_median_filtering(obj, window_size=7, on_references=FALSE)
	expect_equal(unname(only_obs@expr.data[, 1:3]), expr[, 1:3])
	expect_equal(unname(only_obs@expr.data[, 4:6]), matrix(median(expr[, 4:6]), 3, 3))
	only_ref <- infercnv:::apply_median_filtering(obj, window_size=7, on_observations=FALSE)
	expect_equal(unname(only_ref@expr.data[, 1:3]), matrix(5, 3, 3))
	expect_equal(unname(only_ref@expr.data[, 4:6]), expr[, 4:6])
	neither <- infercnv:::apply_median_filtering(obj, window_size=7, on_observations=FALSE, on_references=FALSE)
	expect_equal(neither@expr.data, obj@expr.data)
})

## ---------------------------------------------------------------------------------------------
## masking of non-DE genes (R/inferCNV_mask_non_DE.R)
## ---------------------------------------------------------------------------------------------

## 4 genes x 10 cells; cells 1-5 are the reference (n1), cells 6-10 the only tumor subcluster
## g1, g2: tumor clearly above the reference (no ties, no overlap)
## g3, g4: reference and tumor values interleave, so no DE
opsrun_de_obj <- function() {
	expr <- rbind(g1=c(1, 2, 3, 4, 5,   11, 12, 13, 14, 15),
				  g2=c(2, 3, 4, 5, 6,   21, 22, 23, 24, 25),
				  g3=c(1, 3, 5, 7, 9,   2, 4, 6, 8, 10),
				  g4=c(2, 4, 6, 8, 10,  1, 3, 5, 7, 9))
	obj <- opsrun_obj(expr, ref=list(n1=1:5), obs=list(tum=6:10))
	obj@tumor_subclusters <- list(subclusters=list(tum=list(tum_s1=6:10)))
	obj
}

test_that("get_DE_genes_basic with the t test reports BH-adjusted p-values per subcluster and normal group", {
	obj <- opsrun_de_obj()
	res <- infercnv:::get_DE_genes_basic(obj, p_val_thresh=0.05, test.use="t")
	expect_equal(names(res), "tum_s1,n1")
	r <- res[["tum_s1,n1"]]
	expect_equal(r$tumor_indices, 6:10)
	expect_equal(r$normal, "n1")
	raw <- sapply(1:4, function(i) t.test(obj@expr.data[i, 1:5], obj@expr.data[i, 6:10])$p.value)
	expect_equal(r$pvals, setNames(p.adjust(raw, "BH"), paste0("g", 1:4)))
	expect_equal(r$de_genes, paste0("g", which(p.adjust(raw, "BH") < 0.05)))
	expect_equal(r$de_genes, c("g1", "g2"))
	## the threshold is strict and applied to the adjusted values
	res_strict <- infercnv:::get_DE_genes_basic(obj, p_val_thresh=min(r$pvals) / 10, test.use="t")
	expect_equal(res_strict[[1]]$de_genes, character(0))
})

test_that("get_DE_genes_basic with wilcoxon gives the exact rank-test p-values for clearly separated ranks", {
	set.seed(1)
	obj <- opsrun_de_obj()
	res <- infercnv:::get_DE_genes_basic(obj, p_val_thresh=0.05, test.use="wilcoxon")
	r <- res[["tum_s1,n1"]]
	## the added noise (sd 1e-4) cannot change any rank here, so the p-values are those of the rank test on
	## the unperturbed values
	raw <- sapply(1:4, function(i) wilcox.test(obj@expr.data[i, 1:5], obj@expr.data[i, 6:10])$p.value)
	expect_equal(unname(raw[1:2]), rep(2 / choose(10, 5), 2))
	expect_equal(r$pvals, setNames(p.adjust(raw, "BH"), paste0("g", 1:4)))
	expect_equal(r$de_genes, c("g1", "g2"))
	expect_equal(r$tumor_indices, 6:10)
})

test_that("get_DE_genes_basic runs every subcluster against every normal group", {
	expr <- rbind(g1=c(1, 2, 3,  50, 51, 52,  11, 12, 13,  14, 15, 16),
				  g2=c(5, 6, 7,  8, 9, 10,   20, 21, 22,  60, 61, 62))
	obj <- opsrun_obj(expr, ref=list(nA=1:3, nB=4:6), obs=list(tum=7:12))
	obj@tumor_subclusters <- list(subclusters=list(tum=list(s1=7:9, s2=10:12)))
	res <- infercnv:::get_DE_genes_basic(obj, test.use="t")
	expect_equal(names(res), c("s1,nA", "s1,nB", "s2,nA", "s2,nB"))
	expect_equal(res[["s2,nB"]]$normal, "nB")
	expect_equal(res[["s2,nB"]]$tumor_indices, 10:12)
	expect_equal(res[["s1,nA"]]$tumor_indices, 7:9)
	for (nm in names(res)) {
		expect_equal(names(res[[nm]]$pvals), c("g1", "g2"))
	}
})

test_that("get_DE_genes_basic with the permutation test flags a clearly shifted gene", {
	skip_if_not_installed("coin")
	obj <- opsrun_de_obj()
	res <- infercnv:::get_DE_genes_basic(obj, p_val_thresh=0.05, test.use="perm")
	r <- res[["tum_s1,n1"]]
	expect_true(all(c("g1", "g2") %in% r$de_genes))
	expect_false(any(c("g3", "g4") %in% r$de_genes))
})

test_that(".mask_DE_genes masks genes by the number of normal groups they are DE against", {
	## 12 cells: normals n1=1:2, n2=3:4, n3=5:6 ; tumor subclusters sA = 7:9, sB = 10:12
	expr <- matrix(seq_len(48), nrow=4)
	rownames(expr) <- paste0("g", 1:4)
	obj <- opsrun_obj(expr, ref=list(n1=1:2, n2=3:4, n3=5:6), obs=list(tum=7:12))
	de <- list("sA,n1"=list(tumor_indices=7:9, normal="n1", de_genes=c("g1", "g2", "g3")),
			   "sA,n2"=list(tumor_indices=7:9, normal="n2", de_genes=c("g1", "g2")),
			   "sA,n3"=list(tumor_indices=7:9, normal="n3", de_genes=c("g1")),
			   "sB,n1"=list(tumor_indices=10:12, normal="n1", de_genes=c("g2")),
			   "sB,n2"=list(tumor_indices=10:12, normal="n2", de_genes=c("g2")),
			   "sB,n3"=list(tumor_indices=10:12, normal="n3", de_genes=character(0)))
	## DE counts: sA g1=3 g2=2 g3=1 g4=0 ; sB g2=2 and 0 elsewhere
	any_res <- infercnv:::.mask_DE_genes(obj, de, mask_val=-99, require_DE_all_normals="any", min_cluster_size_mask=3)
	expected_any <- matrix(FALSE, 4, 12)
	expected_any[4, 7:9] <- TRUE
	expected_any[c(1, 3, 4), 10:12] <- TRUE
	expect_equal(unname(any_res@expr.data == -99), expected_any)
	most_res <- infercnv:::.mask_DE_genes(obj, de, mask_val=-99, require_DE_all_normals="most", min_cluster_size_mask=3)
	expected_most <- matrix(FALSE, 4, 12)
	expected_most[c(3, 4), 7:9] <- TRUE
	expected_most[c(1, 3, 4), 10:12] <- TRUE
	expect_equal(unname(most_res@expr.data == -99), expected_most)
	all_res <- infercnv:::.mask_DE_genes(obj, de, mask_val=-99, require_DE_all_normals="all", min_cluster_size_mask=3)
	expected_all <- matrix(FALSE, 4, 12)
	expected_all[c(2, 3, 4), 7:9] <- TRUE
	expected_all[, 10:12] <- TRUE
	expect_equal(unname(all_res@expr.data == -99), expected_all)
	## everything that is not masked keeps its value, reference cells are never touched
	for (res in list(any_res, most_res, all_res)) {
		keep <- res@expr.data != -99
		expect_equal(res@expr.data[keep], expr[keep])
		expect_equal(res@expr.data[, 1:6], obj@expr.data[, 1:6])
	}
})

test_that(".mask_DE_genes never masks subclusters smaller than min_cluster_size_mask", {
	expr <- matrix(as.numeric(seq_len(24)), nrow=4)
	obj <- opsrun_obj(expr, ref=list(n1=1:2), obs=list(tum=3:6))
	de <- list("sA,n1"=list(tumor_indices=3:4, normal="n1", de_genes=character(0)),
			   "sB,n1"=list(tumor_indices=5:6, normal="n1", de_genes=character(0)))
	res <- infercnv:::.mask_DE_genes(obj, de, mask_val=-1, require_DE_all_normals="any")
	expect_equal(res@expr.data, obj@expr.data)
	res2 <- infercnv:::.mask_DE_genes(obj, de, mask_val=-1, require_DE_all_normals="any", min_cluster_size_mask=2)
	expect_equal(unname(res2@expr.data[, 3:6]), matrix(-1, 4, 4))
	expect_equal(res2@expr.data[, 1:2], obj@expr.data[, 1:2])
})

test_that(".mask_DE_genes stops on an unknown require_DE_all_normals", {
	obj <- opsrun_obj(matrix(1:12, nrow=3), ref=list(n1=1:2), obs=list(tum=3:4))
	expect_error(infercnv:::.mask_DE_genes(obj, list(), mask_val=0, require_DE_all_normals="some"),
				 "not recognizing require_DE_all_normals=some")
})

test_that("mask_non_DE_genes_basic sets the non-DE genes of tumor cells to the global mean", {
	set.seed(2)
	obj <- opsrun_de_obj()
	res <- infercnv:::mask_non_DE_genes_basic(obj, p_val_thresh=0.05, test.use="wilcoxon")
	centre <- mean(obj@expr.data)
	expected <- obj@expr.data
	expected[c("g3", "g4"), 6:10] <- centre
	expect_equal(res@expr.data, expected)
	## explicit centre value
	res2 <- infercnv:::mask_non_DE_genes_basic(obj, p_val_thresh=0.05, test.use="t", center_val=0)
	expected2 <- obj@expr.data
	expected2[c("g3", "g4"), 6:10] <- 0
	expect_equal(res2@expr.data, expected2)
})

## ---------------------------------------------------------------------------------------------
## tumor subclustering (R/inferCNV_tumor_subclusters.R)
## ---------------------------------------------------------------------------------------------

## three tight groups of cells with different 8-gene patterns (A, B, C); C is closer to A than to B
opsrun_pattern_matrix <- function(groups) {
	base <- list(A=c(0, 0, 0, 0, 10, 10, 10, 10), B=c(10, 10, 10, 10, 0, 0, 0, 0), C=c(2, 2, 2, 2, 8, 8, 8, 8))
	m <- sapply(seq_along(groups), function(i) base[[groups[i]]] + ((i * 7) %% 5) * 0.01)
	colnames(m) <- paste0("c", seq_along(groups))
	rownames(m) <- paste0("g", seq_len(nrow(m)))
	m
}

## two groups of 12 cells with opposite up/down patterns over `ng` genes, no per-cell offset
opsrun_two_group_matrix <- function(nA=12, nB=12, ng=20) {
	n <- nA + nB
	m <- matrix(0, ng, n)
	for (j in seq_len(n)) {
		a <- j <= nA
		m[, j] <- c(rep(if (a) 3 else -3, ng / 2), rep(if (a) -3 else 3, ng / 2)) + sin(seq_len(ng) * j) * 0.3
	}
	colnames(m) <- paste0("c", seq_len(n))
	rownames(m) <- paste0("g", seq_len(ng))
	m
}

opsrun_named_idx <- function(m) setNames(seq_len(ncol(m)), colnames(m))

test_that(".single_tumor_subclustering with partition 'none' keeps all cells in one subcluster, ordered by the tree", {
	m <- opsrun_pattern_matrix(c("B", "B", "B", "A", "A", "A", "C", "C", "C"))
	idx <- opsrun_named_idx(m)
	res <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=0.1, hclust_method="ward.D2", partition_method="none")
	ref_hc <- stats::hclust(stats::dist(t(m)), method="ward.D2")
	expect_equal(names(res$subclusters), "T_s1")
	expect_equal(res$subclusters[["T_s1"]], idx[ref_hc$order])
	expect_s3_class(res$hc, "hclust")
	expect_equal(res$hc$height, ref_hc$height)
	expect_equal(res$hc$merge, ref_hc$merge)
	expect_equal(res$hc$labels, colnames(m))
})

test_that(".single_tumor_subclustering with 'qnorm' cuts the tree at the normal quantile of its heights", {
	m <- opsrun_pattern_matrix(c("B", "B", "B", "A", "A", "A", "C", "C", "C"))
	idx <- opsrun_named_idx(m)
	hc <- stats::hclust(stats::dist(t(m)), method="ward.D2")
	grps <- stats::cutree(hc, h=qnorm(0.9, mean=mean(hc$height), sd=sd(hc$height)))
	## the cut gives {B} and {A, C}
	expect_equal(unname(grps), c(1, 1, 1, 2, 2, 2, 2, 2, 2))
	res <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=0.1, hclust_method="ward.D2", partition_method="qnorm")
	expect_equal(names(res$subclusters), c("T_s1", "T_s2"))
	expect_equal(sort(names(res$subclusters$T_s1)), paste0("c", 1:3))
	expect_equal(sort(names(res$subclusters$T_s2)), paste0("c", 4:9))
	## every cell once, index and name agree
	all_idx <- unlist(res$subclusters, use.names=FALSE)
	expect_equal(sort(all_idx), 1:9)
	for (s in res$subclusters) expect_equal(names(s), colnames(m)[s])
	## within a subcluster the cells follow the tree order
	expect_equal(unname(unlist(res$subclusters)), unname(idx[hc$order]))
})

test_that(".single_tumor_subclustering with 'pheight' cuts at p_val times the maximum height", {
	## cell 1-4 pattern A, 5-8 pattern B
	m <- opsrun_pattern_matrix(c("A", "A", "A", "A", "B", "B", "B", "B"))
	idx <- opsrun_named_idx(m)
	hc <- stats::hclust(stats::dist(t(m)), method="ward.D2")
	res <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=0.5, hclust_method="ward.D2", partition_method="pheight")
	expect_equal(names(res$subclusters), c("T_s1", "T_s2"))
	expect_equal(sort(unname(res$subclusters$T_s1)), 1:4)
	expect_equal(sort(unname(res$subclusters$T_s2)), 5:8)
	## a cut above the maximum height gives one cluster
	res1 <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=1.5, hclust_method="ward.D2", partition_method="pheight")
	expect_equal(names(res1$subclusters), "T_s1")
	expect_equal(sort(unname(res1$subclusters$T_s1)), 1:8)
})

test_that(".single_tumor_subclustering with 'qgamma' cuts at the gamma quantile of the heights", {
	skip_if_not_installed("fitdistrplus")
	m <- opsrun_pattern_matrix(c("B", "B", "B", "A", "A", "A", "C", "C", "C"))
	idx <- opsrun_named_idx(m)
	hc <- stats::hclust(stats::dist(t(m)), method="ward.D2")
	fit <- fitdistrplus::fitdist(hc$height, "gamma")
	cut <- qgamma(0.9, shape=fit$estimate[[1]], rate=fit$estimate[[2]])
	expected <- stats::cutree(hc, h=cut)
	res <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=0.1, hclust_method="ward.D2", partition_method="qgamma")
	expect_equal(length(res$subclusters), length(unique(expected)))
	expect_equal(sort(unname(sapply(res$subclusters, length))), sort(as.vector(table(expected))))
})

test_that(".single_tumor_subclustering with two or fewer cells returns a single subcluster without a tree", {
	m <- opsrun_pattern_matrix(c("A", "B"))
	idx <- opsrun_named_idx(m)
	res <- infercnv:::.single_tumor_subclustering("T", idx, m, p_val=0.1, hclust_method="ward.D2", partition_method="qnorm")
	expect_null(res$hc)
	expect_equal(res$subclusters, list(T_s1=idx))
})

test_that(".find_DE_stat_significance returns the Welch t-test p-value of every gene", {
	normal <- rbind(g1=c(1, 2, 3, 4), g2=c(5, 5, 5, 6), g3=c(1, 1, 1, 2))
	tumor <- rbind(g1=c(10, 11, 12, 13), g2=c(5, 6, 5, 5), g3=c(2, 2, 2, 3))
	res <- infercnv:::.find_DE_stat_significance(normal, tumor)
	expect_equal(res, sapply(1:3, function(i) t.test(normal[i, ], tumor[i, ])$p.value))
	expect_lt(res[1], 0.001)
	expect_gt(res[2], 0.5)
	## constant vs constant: the test fails and NA is reported
	expect_true(is.na(infercnv:::.find_DE_stat_significance(rbind(c(1, 1, 1)), rbind(c(1, 1, 1)))))
})

opsrun_subcluster_obj <- function() {
	## 20 genes on 2 chromosomes ; cells 1-4 reference, 5-12 tumor 'tA', 13-20 tumor 'tB'
	set.seed(11)
	m <- matrix(0, 20, 20)
	for (j in 1:20) m[, j] <- 2 + sin(seq_len(20) * j) * 0.2
	m[1:10, 5:12] <- m[1:10, 5:12] + 3
	m[11:20, 13:20] <- m[11:20, 13:20] - 1
	colnames(m) <- paste0("c", 1:20)
	rownames(m) <- paste0("g", 1:20)
	opsrun_obj(m, chr=rep(c("chr1", "chr2"), each=10), ref=list(r=1:4), obs=list(tA=5:12, tB=13:20))
}

test_that("define_signif_tumor_subclusters with partition 'none' makes one subcluster and one tree per group", {
	obj <- opsrun_subcluster_obj()
	res <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.1, partition_method="none", z_score_filter=0,
													  cluster_by_groups=TRUE)
	expect_equal(length(res), 2)
	expect_null(res[[2]])
	o <- res[[1]]
	expect_equal(names(o@tumor_subclusters$subclusters), c("tA", "tB", "r"))
	expect_equal(names(o@tumor_subclusters$hc), c("tA", "tB", "r"))
	for (g in c("tA", "tB", "r")) {
		cells <- c(tA=5:12, tB=13:20, r=1:4)[grepl(paste0("^", g), names(c(tA=5:12, tB=13:20, r=1:4)))]
		cells <- unname(cells)
		sc <- o@tumor_subclusters$subclusters[[g]]
		expect_equal(names(sc), paste0(g, "_s1"))
		expect_equal(sort(unname(sc[[1]])), sort(cells))
		expect_equal(names(sc[[1]]), colnames(obj@expr.data)[sc[[1]]])
		ref_hc <- stats::hclust(stats::dist(t(obj@expr.data[, cells])), method="ward.D2")
		expect_equal(o@tumor_subclusters$hc[[g]]$height, ref_hc$height)
		expect_equal(unname(sc[[1]]), cells[ref_hc$order])
	}
	## data and layout stay as they were
	expect_identical(o@expr.data, obj@expr.data)
	expect_identical(o@observation_grouped_cell_indices, obj@observation_grouped_cell_indices)
})

test_that("define_signif_tumor_subclusters with cluster_by_groups=FALSE merges all observations", {
	obj <- opsrun_subcluster_obj()
	res <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.1, partition_method="none", z_score_filter=0,
													  cluster_by_groups=FALSE)
	o <- res[[1]]
	expect_equal(names(o@tumor_subclusters$subclusters), c("all_observations", "r"))
	expect_equal(sort(unname(o@tumor_subclusters$subclusters$all_observations[[1]])), 5:20)
	expect_equal(names(o@tumor_subclusters$subclusters$all_observations), "all_observations_s1")
})

test_that("define_signif_tumor_subclusters with pheight separates tumor cells with different profiles", {
	obj <- opsrun_subcluster_obj()
	## merge the two tumor groups: tA (up on chr1) and tB (down on chr2) are well separated
	obj@observation_grouped_cell_indices <- list(t=5:20)
	res <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.5, partition_method="pheight", z_score_filter=0,
													  cluster_by_groups=TRUE)
	sc <- res[[1]]@tumor_subclusters$subclusters$t
	expect_equal(length(sc), 2)
	membership <- sapply(sc, function(s) paste(sort(unname(s)), collapse=","))
	expect_equal(sort(unname(membership)), sort(c(paste(5:12, collapse=","), paste(13:20, collapse=","))))
})

test_that("define_signif_tumor_subclusters drops reference z-score outlier genes before clustering", {
	## 6 genes, reference cells 1-4 ; gene g1 is an outlier (mean |z| over the reference >= 0.8)
	m <- rbind(g1=c(100, 100, 100, 100, 0, 50, 100, 20),
			   g2=c(1, 2, 3, 4, 3, 2, 2, 4),
			   g3=c(4, 3, 2, 1, 1, 4, 3, 2),
			   g4=c(2, 1, 4, 3, 2, 3, 4, 1),
			   g5=c(3, 4, 1, 2, 4, 1, 1, 3),
			   g6=c(1, 2, 3, 4, 2, 2, 4, 3))
	colnames(m) <- paste0("c", 1:8)
	obj <- opsrun_obj(m, ref=list(r=1:4), obs=list(t=5:8))
	ref_matrix <- m[, 1:4]
	z <- (ref_matrix - mean(ref_matrix)) / sd(ref_matrix)
	outliers <- which(apply(abs(z), 1, mean) >= 0.8)
	expect_equal(unname(outliers), 1L)
	res <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.1, partition_method="none", z_score_filter=0.8)
	ref_hc <- stats::hclust(stats::dist(t(m[-1, 5:8])), method="ward.D2")
	expect_equal(res[[1]]@tumor_subclusters$hc$t$height, ref_hc$height)
	## without the filter the outlier gene dominates the distances
	res0 <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.1, partition_method="none", z_score_filter=0)
	ref_hc0 <- stats::hclust(stats::dist(t(m[, 5:8])), method="ward.D2")
	expect_equal(res0[[1]]@tumor_subclusters$hc$t$height, ref_hc0$height)
	expect_false(isTRUE(all.equal(ref_hc$height, ref_hc0$height)))
})

test_that("define_signif_tumor_subclusters with restrict_to_DE_genes clusters on genes that differ from the reference", {
	obj <- opsrun_subcluster_obj()
	res <- infercnv:::define_signif_tumor_subclusters(obj, p_val=0.01, partition_method="none", z_score_filter=0,
													  restrict_to_DE_genes=TRUE)
	o <- res[[1]]
	## the tumor group tA differs from the reference on chr1 only; its tree is built on the DE genes
	ref_cells <- 1:4
	tA <- 5:12
	p <- sapply(1:20, function(i) t.test(obj@expr.data[i, ref_cells], obj@expr.data[i, tA])$p.value)
	de <- which(p < 0.01)
	expect_true(length(de) > 0)
	ref_hc <- stats::hclust(stats::dist(t(obj@expr.data[de, tA, drop=FALSE])), method="ward.D2")
	expect_equal(o@tumor_subclusters$hc$tA$height, ref_hc$height)
})

test_that(".single_tumor_leiden_subclustering keeps tiny groups as one subcluster", {
	m <- opsrun_two_group_matrix()
	idx <- opsrun_named_idx(m)
	res2 <- infercnv:::.single_tumor_leiden_subclustering("T", idx[1:2], m[, 1:2], k_nn=5, leiden_resolution=1,
														  leiden_method="simple", leiden_function="modularity",
														  hclust_method="ward.D2")
	expect_null(res2$hc)
	expect_equal(res2$subclusters, list(T_s1=idx[1:2]))
	## fewer cells than k_nn: a single subcluster named after the group, with a tree
	res5 <- infercnv:::.single_tumor_leiden_subclustering("T", idx[1:5], m[, 1:5], k_nn=5, leiden_resolution=1,
														  leiden_method="simple", leiden_function="modularity",
														  hclust_method="ward.D2")
	expect_equal(names(res5$subclusters), "T")
	expect_equal(res5$subclusters$T, idx[1:5])
	expect_s3_class(res5$hc, "hclust")
	expect_equal(res5$hc$height, stats::hclust(stats::dist(t(m[, 1:5])), "ward.D2")$height)
})

test_that(".single_tumor_leiden_subclustering never mixes two separated groups of cells", {
	set.seed(3)
	m <- opsrun_two_group_matrix()
	idx <- opsrun_named_idx(m)
	res <- infercnv:::.single_tumor_leiden_subclustering("T", idx, m, k_nn=5, leiden_resolution=1,
														 leiden_method="simple", leiden_function="modularity",
														 hclust_method="ward.D2")
	expect_true(length(res$subclusters) >= 2)
	expect_equal(sort(unname(unlist(res$subclusters))), 1:24)
	for (s in res$subclusters) {
		expect_true(all(s <= 12) || all(s >= 13))
		expect_equal(names(s), colnames(m)[s])
	}
	expect_true(all(grepl("^T_s[0-9]+$", names(res$subclusters))))
	## the tree has all cells as leaves
	expect_s3_class(res$hc, "hclust")
	expect_equal(sort(res$hc$labels), sort(colnames(m)))
})

test_that(".single_tumor_leiden_subclustering with the PCA preprocessing never mixes two separated groups", {
	skip_if_not_installed("Seurat")
	skip_if(utils::packageVersion("Seurat") < "5.0.0")
	set.seed(3)
	m <- opsrun_two_group_matrix()
	idx <- opsrun_named_idx(m)
	invisible(utils::capture.output(suppressMessages(suppressWarnings(
		out <- infercnv:::.single_tumor_leiden_subclustering("T", idx, m, k_nn=5, leiden_resolution=1,
															 leiden_method="PCA", leiden_function="modularity",
															 hclust_method="ward.D2")))))
	expect_equal(sort(unname(unlist(out$subclusters))), 1:24)
	for (s in out$subclusters) expect_true(all(s <= 12) || all(s >= 13))
})

test_that(".leiden_simple_snn returns one label per cell and respects disconnected groups", {
	set.seed(4)
	m <- opsrun_two_group_matrix()
	for (fn in c("modularity", "CPM")) {
		p <- infercnv:::.leiden_simple_snn(m, k_nn=5, resolution_parameter=1, objective_function=fn)
		expect_equal(length(p), 24)
		## k_nn = 5 nearest neighbours of a cell are always in its own group
		for (k in unique(p)) {
			in_k <- which(p == k)
			expect_true(all(in_k <= 12) || all(in_k >= 13))
		}
	}
})

test_that(".whole_dataset_leiden_subclustering_per_chr keeps small groups and single-cluster cases as they are", {
	m <- opsrun_two_group_matrix()
	chrs <- factor(rep(c("chr1", "chr2"), each=10), levels=c("chr1", "chr2", "chr3"))
	groups <- list(a=setNames(1:12, colnames(m)[1:12]), b=setNames(13:24, colnames(m)[13:24]),
				   small=setNames(c(1L, 13L), colnames(m)[c(1, 13)]))
	## k_nn larger than every group: each chromosome keeps the groups whole
	res <- infercnv:::.whole_dataset_leiden_subclustering_per_chr(m, groups, chrs, k_nn=15, leiden_resolution=1,
																  leiden_method="simple", leiden_function="modularity")
	expect_equal(names(res), c("chr1", "chr2", "chr3"))
	for (ch in c("chr1", "chr2")) {
		expect_equal(res[[ch]], groups)
	}
})

# pins the current output
test_that(".whole_dataset_leiden_subclustering_per_chr gives a chromosome level without genes all cells for every group", {
	m <- opsrun_two_group_matrix()
	chrs <- factor(rep(c("chr1", "chr2"), each=10), levels=c("chr1", "chr2", "chr3"))
	groups <- list(a=setNames(1:12, colnames(m)[1:12]), b=setNames(13:24, colnames(m)[13:24]))
	res <- infercnv:::.whole_dataset_leiden_subclustering_per_chr(m, groups, chrs, k_nn=15, leiden_resolution=1,
																  leiden_method="simple", leiden_function="modularity")
	expect_equal(names(res$chr3), c("a", "b"))
	for (g in names(res$chr3)) expect_equal(unname(res$chr3[[g]]), 1:24)
})

test_that(".whole_dataset_leiden_subclustering_per_chr clusters each group of each chromosome", {
	set.seed(5)
	m <- opsrun_two_group_matrix()
	chrs <- factor(rep(c("chr1", "chr2"), each=10), levels=c("chr1", "chr2"))
	groups <- list(a=setNames(1:12, colnames(m)[1:12]), b=setNames(13:24, colnames(m)[13:24]))
	res <- infercnv:::.whole_dataset_leiden_subclustering_per_chr(m, groups, chrs, k_nn=5, leiden_resolution=1,
																  leiden_method="simple", leiden_function="modularity")
	for (ch in c("chr1", "chr2")) {
		sub <- res[[ch]]
		expect_true(all(grepl("^[ab]_s[0-9]+$", names(sub))))
		a_cells <- unlist(sub[grepl("^a_", names(sub))], use.names=FALSE)
		b_cells <- unlist(sub[grepl("^b_", names(sub))], use.names=FALSE)
		expect_equal(sort(a_cells), 1:12)
		expect_equal(sort(b_cells), 13:24)
		for (s in sub) expect_equal(names(s), colnames(m)[s])
	}
})

test_that("add_single_branch_to_phylo adds one tip at the root of a rooted ultrametric tree", {
	m <- opsrun_two_group_matrix()[, c(1, 2, 13)]
	tr <- ape::as.phylo(stats::hclust(stats::dist(t(m)), "complete"))
	nt <- infercnv:::add_single_branch_to_phylo(tr, "new")
	expect_equal(nt$tip.label, c("new", tr$tip.label))
	expect_equal(ape::Ntip(nt), 4)
	expect_equal(nt$Nnode, tr$Nnode + 1)
	expect_true(ape::is.rooted(nt))
	h <- max(ape::node.depth.edgelength(tr)[seq_len(ape::Ntip(tr))])
	## new root splits into the new tip (length h + 1) and the old tree (root edge 1)
	expect_equal(nt$edge.length[1:2], c(h + 1, 1))
	expect_equal(nt$edge.length[-(1:2)], tr$edge.length)
	## distances between the old tips are unchanged, the new tip is 2h + 2 from every old tip
	d_old <- ape::cophenetic.phylo(tr)
	d_new <- ape::cophenetic.phylo(nt)
	expect_equal(d_new[rownames(d_old), colnames(d_old)], d_old)
	expect_equal(unname(d_new["new", rownames(d_old)]), rep(2 * h + 2, 3))
})

## ---------------------------------------------------------------------------------------------
## random-tree subclustering
## ---------------------------------------------------------------------------------------------

test_that(".get_tree_height_via_ecdf returns the quantile of the random maximum heights", {
	params <- list(ecdf=ecdf(1:10))
	expect_equal(unname(infercnv:::.get_tree_height_via_ecdf(0.1, params)), 9.1)
	expect_equal(unname(infercnv:::.get_tree_height_via_ecdf(0.5, params)), 5.5)
})

test_that(".parameterize_random_cluster_heights reports the observed tree and 100 permuted maximum heights", {
	set.seed(6)
	m <- opsrun_two_group_matrix(ng=30)
	res <- suppressMessages(infercnv:::.parameterize_random_cluster_heights(m, "ward.D2", plot=FALSE))
	ref_hc <- stats::hclust(stats::dist(t(m)), "ward.D2")
	expect_equal(res$h_obs$height, ref_hc$height)
	expect_equal(res$max_h, max(ref_hc$height))
	expect_equal(length(res$rand_max_height_dist), 100)
	## permuting every gene across the cells destroys the A/B structure, so the observed maximum is far above
	expect_true(all(res$rand_max_height_dist < res$max_h))
	expect_equal(res$ecdf(max(res$rand_max_height_dist)), 1)
	expect_equal(res$ecdf(res$max_h), 1)
	expect_s3_class(res$h_rand_ex, "hclust")
})

test_that(".parameterize_random_cluster_heights_smoothed_trees uses smoothed, median-centred data", {
	set.seed(7)
	m <- opsrun_two_group_matrix(ng=30)
	res <- suppressWarnings(suppressMessages(
		infercnv:::.parameterize_random_cluster_heights_smoothed_trees(m, "ward.D2", window_size=5)))
	sm <- apply(m, 2, caTools::runmean, k=5)
	sm <- sweep(sm, 2, apply(sm, 2, median), "-")
	ref_hc <- stats::hclust(stats::dist(t(sm)), "ward.D2")
	expect_equal(res$h_obs$height, ref_hc$height)
	expect_equal(res$max_h, max(ref_hc$height))
	expect_equal(length(res$rand_max_height_dist), 100)
	expect_true(all(res$rand_max_height_dist < res$max_h))
	expect_equal(res$ecdf(res$max_h), 1)
})

test_that(".partition_by_random_smoothed_trees splits two well separated groups at the first level", {
	set.seed(8)
	m <- opsrun_two_group_matrix(ng=30)
	invisible(utils::capture.output(suppressWarnings(suppressMessages(
		grps <- infercnv:::.partition_by_random_smoothed_trees("T", m, "ward.D2", p_val=0.1, window_size=5,
															   max_recursion_depth=1, min_cluster_size_recurse=10)))))
	expect_equal(names(grps), colnames(m))
	expect_equal(unname(grps), rep(c("T.1.1", "T.1.2"), each=12))
})

test_that(".partition_by_random_smoothed_trees stops at the maximum recursion depth", {
	set.seed(8)
	m <- opsrun_two_group_matrix(ng=30)
	## depth 0 is below the first level, so nothing is evaluated and all cells stay in the root clade
	invisible(utils::capture.output(suppressWarnings(suppressMessages(
		grps <- infercnv:::.single_tumor_subclustering_recursive_random_smoothed_trees(
			m, "ward.D2", 0.1, setNames(rep("T.1", 24), colnames(m)), window_size=5,
			max_recursion_depth=0, min_cluster_size_recurse=10)))))
	expect_equal(unname(grps), rep("T.1", 24))
})

test_that("define_signif_tumor_subclusters_via_random_smooothed_trees splits a tumor group on its profile", {
	set.seed(9)
	ng <- 30
	## log-scale data: 8 reference cells near 0, a tumor with two groups of 12 cells with opposite patterns
	m <- matrix(0, ng, 32)
	for (j in 1:32) m[, j] <- sin(seq_len(ng) * j) * 0.1
	pat <- c(rep(3, ng / 2), rep(-3, ng / 2))
	m[, 9:20] <- m[, 9:20] + pat
	m[, 21:32] <- m[, 21:32] - pat
	colnames(m) <- paste0("c", 1:32)
	rownames(m) <- paste0("g", 1:ng)
	obj <- opsrun_obj(m, ref=list(r=1:8), obs=list(t=9:32))
	invisible(utils::capture.output(suppressWarnings(suppressMessages(
		res <- infercnv:::define_signif_tumor_subclusters_via_random_smooothed_trees(
			obj, p_val=0.1, hclust_method="ward.D2", cluster_by_groups=TRUE, window_size=5,
			max_recursion_depth=1, min_cluster_size_recurse=10)))))
	sc <- res@tumor_subclusters$subclusters
	expect_equal(names(sc), c("t", "r"))
	expect_equal(names(sc$t), c("t.1.1", "t.1.2"))
	expect_equal(unname(sc$t[["t.1.1"]]), 9:20)
	expect_equal(unname(sc$t[["t.1.2"]]), 21:32)
	expect_equal(names(sc$t[["t.1.1"]]), colnames(m)[9:20])
	## the reference group is partitioned too: every reference cell is in exactly one subcluster
	expect_true(all(grepl("^r\\.", names(sc$r))))
	expect_equal(sort(unname(unlist(sc$r))), 1:8)
	expect_equal(names(res@tumor_subclusters$hc), c("t", "r"))
	expect_s3_class(res@tumor_subclusters$hc$t, "hclust")
	expect_equal(sort(res@tumor_subclusters$hc$t$labels), sort(colnames(m)[9:32]))
	## the expression data of the object is not altered
	expect_identical(res@expr.data, obj@expr.data)
})

## ---------------------------------------------------------------------------------------------
## run(): early checks, saved step objects and a tiny end-to-end call
## ---------------------------------------------------------------------------------------------

## 200 genes (100 on chr1, 100 on chr2) x 30 cells from the example data: 10 reference cells
## (Oligodendrocytes) and 10 cells of each of two malignant samples
opsrun_tiny_cache <- new.env()

opsrun_make_tiny_obj <- function() {
	if (!is.null(opsrun_tiny_cache$obj)) return(opsrun_tiny_cache$obj)
	ext <- function(f) system.file("extdata", f, package="infercnv")
	counts <- read.table(gzfile(ext("oligodendroglioma_expression_downsampled.counts.matrix.gz")),
						 header=TRUE, sep="\t", row.names=1, check.names=FALSE)
	annot <- read.table(ext("oligodendroglioma_annotations_downsampled.txt"), sep="\t", row.names=1,
						stringsAsFactors=FALSE)
	genes <- read.table(ext("gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt"), sep="\t", row.names=1,
						stringsAsFactors=FALSE)
	pick <- function(label) rownames(annot)[annot[, 1] == label][1:10]
	cells <- c(pick("Oligodendrocytes (non-malignant)"), pick("malignant_93"), pick("malignant_97"))
	gg <- genes[genes[, 1] %in% c("chr1", "chr2") & rownames(genes) %in% rownames(counts), ]
	gg <- gg[order(gg[, 1], gg[, 2]), ]
	sel <- c(rownames(gg)[gg[, 1] == "chr1"][round(seq(1, 1000, length.out=100))],
			 rownames(gg)[gg[, 1] == "chr2"][round(seq(1, 700, length.out=100))])
	td <- tempfile("opsrun_tiny_")
	dir.create(td)
	write.table(counts[sel, cells], file.path(td, "counts.txt"), sep="\t", quote=FALSE)
	write.table(annot[cells, , drop=FALSE], file.path(td, "annot.txt"), sep="\t", quote=FALSE, col.names=FALSE)
	write.table(genes[sel, ], file.path(td, "genes.txt"), sep="\t", quote=FALSE, col.names=FALSE)
	obj <- CreateInfercnvObject(file.path(td, "counts.txt"), file.path(td, "genes.txt"), file.path(td, "annot.txt"),
								ref_group_names="Oligodendrocytes (non-malignant)",
								min_max_counts_per_cell=c(0, Inf))
	unlink(td, recursive=TRUE)
	opsrun_tiny_cache$obj <- obj
	obj
}

## genes passing run()'s two filters: mean count >= cutoff and at least min_cells cells above zero
opsrun_expected_genes <- function(counts, cutoff, min_cells) {
	keep <- rowMeans(counts) >= cutoff
	counts <- counts[keep, , drop=FALSE]
	rownames(counts)[rowSums(counts > 0) >= min_cells]
}

test_that("run() stops early on invalid arguments", {
	obj <- opsrun_make_tiny_obj()
	expect_error(infercnv::run(obj, out_dir=NULL, HMM=FALSE, no_plot=TRUE), "out_dir is NULL")
	expect_error(infercnv::run(obj, out_dir=file.path(tempdir(), "opsrun_bad"), HMM=TRUE, HMM_type="i6",
							   smooth_method="coordinates", no_plot=TRUE),
				 "Incompatible HMM mode and smoothing method")
	expect_false(dir.exists(file.path(tempdir(), "opsrun_bad")))
})

test_that("run() with up_to_step returns early and saves only the steps it ran", {
	obj <- opsrun_make_tiny_obj()
	counts <- obj@count.data
	od0 <- file.path(tempfile("opsrun_up0_"))
	r0 <- infercnv::run(obj, cutoff=1, out_dir=od0, HMM=FALSE, no_plot=TRUE, num_threads=1, up_to_step=0)
	expect_true(dir.exists(od0))
	expect_equal(list.files(od0), character(0))
	expect_equal(dim(r0@expr.data), dim(obj@expr.data))
	od2 <- file.path(tempfile("opsrun_up2_"))
	r2 <- infercnv::run(obj, cutoff=1, out_dir=od2, HMM=FALSE, no_plot=TRUE, num_threads=1, up_to_step=2)
	expect_equal(sort(list.files(od2)), c("01_incoming_data.infercnv_obj", "02_reduced_by_cutoff.infercnv_obj"))
	genes <- opsrun_expected_genes(counts, cutoff=1, min_cells=3)
	expect_equal(rownames(r2@expr.data), genes)
	expect_equal(dim(r2@expr.data), c(length(genes), 30))
	expect_equal(rownames(r2@gene_order), genes)
	expect_equal(r2@expr.data, counts[genes, ])
	step1 <- readRDS(file.path(od2, "01_incoming_data.infercnv_obj"))
	expect_equal(dim(step1@expr.data), c(200, 30))
	step2 <- readRDS(file.path(od2, "02_reduced_by_cutoff.infercnv_obj"))
	expect_equal(rownames(step2@expr.data), genes)
	## step 3 normalizes by depth: every cell sums to the median library size of the filtered matrix
	od3 <- file.path(tempfile("opsrun_up3_"))
	r3 <- infercnv::run(obj, cutoff=1, out_dir=od3, HMM=FALSE, no_plot=TRUE, num_threads=1, up_to_step=3)
	expect_equal(unname(colSums(r3@expr.data)), rep(median(colSums(counts[genes, ])), 30))
	expect_true(file.exists(file.path(od3, "03_normalized_by_depth.infercnv_obj")))
	## save_rds = FALSE writes nothing
	odn <- file.path(tempfile("opsrun_nosave_"))
	infercnv::run(obj, cutoff=1, out_dir=odn, HMM=FALSE, no_plot=TRUE, num_threads=1, up_to_step=4, save_rds=FALSE)
	expect_equal(list.files(odn), character(0))
})

test_that("run() on a tiny example writes every step object and returns the final object", {
	obj <- opsrun_make_tiny_obj()
	counts <- obj@count.data
	od <- file.path(tempfile("opsrun_e2e_"))
	res <- infercnv::run(obj, cutoff=1, out_dir=od, cluster_by_groups=TRUE, denoise=TRUE, HMM=FALSE,
						 no_plot=TRUE, num_threads=1, window_length=21,
						 tumor_subcluster_partition_method="qnorm", tumor_subcluster_pval=0.1,
						 inspect_subclusters=FALSE)
	genes <- opsrun_expected_genes(counts, cutoff=1, min_cells=3)
	expect_s4_class(res, "infercnv")
	expect_equal(dim(res@expr.data), c(length(genes), 30))
	expect_equal(rownames(res@expr.data), genes)
	expect_equal(colnames(res@expr.data), colnames(counts))
	expect_equal(rownames(res@gene_order), genes)

	expected_files <- c("01_incoming_data.infercnv_obj",
						"02_reduced_by_cutoff.infercnv_obj",
						"03_normalized_by_depth.infercnv_obj",
						"04_logtransformed.infercnv_obj",
						"08_remove_ref_avg_from_obs_logFC.infercnv_obj",
						"09_apply_max_centered_expr_threshold.infercnv_obj",
						"10_smoothed_by_chr.infercnv_obj",
						"11_recentered_cells_by_chr.infercnv_obj",
						"12_remove_ref_avg_from_obs_adjust.infercnv_obj",
						"14_invert_log_transform.infercnv_obj",
						"15_tumor_subclusters.qnorm.infercnv_obj",
						"22_denoise.qnorm.NF_NA.SD_1.5.NL_FALSE.infercnv_obj",
						"preliminary.infercnv_obj",
						"run.final.infercnv_obj")
	expect_equal(sort(list.files(od)), sort(expected_files))

	## the final object on disk is the returned object, and step 22 is the same object
	final <- readRDS(file.path(od, "run.final.infercnv_obj"))
	expect_equal(final@expr.data, res@expr.data)
	expect_equal(final@tumor_subclusters, res@tumor_subclusters)
	step22 <- readRDS(file.path(od, "22_denoise.qnorm.NF_NA.SD_1.5.NL_FALSE.infercnv_obj"))
	expect_equal(step22@expr.data, res@expr.data)
	prelim <- readRDS(file.path(od, "preliminary.infercnv_obj"))
	step15 <- readRDS(file.path(od, "15_tumor_subclusters.qnorm.infercnv_obj"))
	expect_equal(prelim@expr.data, step15@expr.data)

	## dimensions of the intermediate objects
	for (f in setdiff(expected_files, "01_incoming_data.infercnv_obj")) {
		o <- readRDS(file.path(od, f))
		expect_equal(dim(o@expr.data), c(length(genes), 30), info=f)
	}
	expect_equal(dim(readRDS(file.path(od, "01_incoming_data.infercnv_obj"))@expr.data), c(200, 30))

	## the run options are recorded on the object
	expect_true(res@options$denoise)
	expect_equal(res@options$window_length, 21)
	expect_equal(res@options$tumor_subcluster_partition_method, "qnorm")
	expect_false(res@options$HMM)

	## subclusters cover every cell of each group exactly once
	sc <- res@tumor_subclusters$subclusters
	expect_equal(sort(names(sc)), sort(c("malignant_93", "malignant_97", "Oligodendrocytes (non-malignant)")))
	for (g in names(sc)) {
		cells <- unlist(c(res@observation_grouped_cell_indices, res@reference_grouped_cell_indices)[g], use.names=FALSE)
		expect_equal(sort(unname(unlist(sc[[g]]))), sort(cells), info=g)
	}
	expect_s3_class(res@tumor_subclusters$hc[["malignant_93"]], "hclust")

	## no plot was requested
	expect_false(any(grepl("\\.png$|\\.pdf$", list.files(od))))
})

test_that("the final expression data of run() equals the step functions applied in run()'s order", {
	obj <- opsrun_make_tiny_obj()
	od <- file.path(tempfile("opsrun_chain_"))
	res <- infercnv::run(obj, cutoff=1, out_dir=od, denoise=TRUE, HMM=FALSE, no_plot=TRUE, num_threads=1,
						 window_length=21, tumor_subcluster_partition_method="qnorm", inspect_subclusters=FALSE,
						 save_rds=FALSE)
	o <- infercnv:::require_above_min_mean_expr_cutoff(obj, 1)
	o <- infercnv:::require_above_min_cells_ref(o, 3)
	o <- infercnv:::normalize_counts_by_seq_depth(o)
	o <- infercnv:::log2xplus1(o)
	o <- infercnv:::subtract_ref_expr_from_obs(o, inv_log=FALSE, use_bounds=TRUE)
	o <- infercnv:::apply_max_threshold_bounds(o, threshold=3)
	o <- infercnv:::smooth_by_chromosome(o, window_length=21, smooth_ends=TRUE)
	o <- infercnv:::center_cell_expr_across_chromosome(o, method="median")
	o <- infercnv:::subtract_ref_expr_from_obs(o, inv_log=FALSE, use_bounds=TRUE)
	o <- infercnv:::invert_log2(o)
	o <- infercnv:::clear_noise_via_ref_mean_sd(o, sd_amplifier=1.5, noise_logistic=FALSE)
	expect_equal(res@expr.data, o@expr.data)
	## no files with save_rds = FALSE except the final one
	expect_equal(list.files(od), "run.final.infercnv_obj")
})

test_that("run() with save_final_rds = FALSE does not write the final object", {
	obj <- opsrun_make_tiny_obj()
	od <- file.path(tempfile("opsrun_nofinal_"))
	infercnv::run(obj, cutoff=1, out_dir=od, HMM=FALSE, no_plot=TRUE, num_threads=1, window_length=21,
				  tumor_subcluster_partition_method="pheight", inspect_subclusters=FALSE, save_rds=FALSE,
				  save_final_rds=FALSE)
	expect_equal(list.files(od), character(0))
})

test_that("run() in resume mode reuses the saved results when the settings are unchanged", {
	obj <- opsrun_make_tiny_obj()
	od <- file.path(tempfile("opsrun_resume_"))
	args <- list(cutoff=1, out_dir=od, denoise=TRUE, HMM=FALSE, no_plot=TRUE, num_threads=1, window_length=21,
				 tumor_subcluster_partition_method="qnorm", inspect_subclusters=FALSE)
	first <- do.call(infercnv::run, c(list(obj), args))
	files <- list.files(od, full.names=TRUE)
	mt <- file.mtime(files)
	Sys.sleep(1.1)
	second <- do.call(infercnv::run, c(list(obj), args))
	expect_equal(second@expr.data, first@expr.data)
	expect_equal(second@tumor_subclusters, first@tumor_subclusters)
	## the numbered step objects were not rewritten
	steps <- grepl("/[0-9][0-9]_", files)
	expect_equal(file.mtime(files)[steps], mt[steps])
})

test_that("run() with several options writes the matching step objects and applies them", {
	obj <- opsrun_make_tiny_obj()
	counts <- obj@count.data
	od <- file.path(tempfile("opsrun_opts_"))
	set.seed(21)
	res <- infercnv::run(obj, cutoff=1, out_dir=od, HMM=FALSE, no_plot=TRUE, num_threads=1,
						 analysis_mode="samples", smooth_method="runmeans", window_length=11,
						 remove_genes_at_chr_ends=TRUE, prune_outliers=TRUE, mask_nonDE_genes=TRUE,
						 test.use="t", denoise=TRUE, noise_filter=0.05, save_final_rds=TRUE)
	expected_files <- c("01_incoming_data.infercnv_obj",
						"02_reduced_by_cutoff.infercnv_obj",
						"03_normalized_by_depth.infercnv_obj",
						"04_logtransformed.infercnv_obj",
						"08_remove_ref_avg_from_obs_logFC.infercnv_obj",
						"09_apply_max_centered_expr_threshold.infercnv_obj",
						"10_smoothed_by_chr.infercnv_obj",
						"11_recentered_cells_by_chr.infercnv_obj",
						"12_remove_ref_avg_from_obs_adjust.infercnv_obj",
						"13_remove_gene_at_chr_ends.infercnv_obj",
						"14_invert_log_transform.infercnv_obj",
						"15_no_subclustering.infercnv_obj",
						"16_remove_outlier.infercnv_obj",
						"21_mask_nonDE.infercnv_obj",
						"22_denoise.NF_0.05.SD_1.5.NL_FALSE.infercnv_obj",
						"preliminary.infercnv_obj",
						"run.final.infercnv_obj")
	expect_equal(sort(list.files(od)), sort(expected_files))
	genes <- opsrun_expected_genes(counts, cutoff=1, min_cells=3)
	chr <- as.character(obj@gene_order[genes, "chr"])
	## window 11 -> 5 genes are removed at both ends of each chromosome
	expect_equal(dim(res@expr.data), c(length(genes) - 10 * length(unique(chr)), 30))
	step12 <- readRDS(file.path(od, "12_remove_ref_avg_from_obs_adjust.infercnv_obj"))
	step13 <- readRDS(file.path(od, "13_remove_gene_at_chr_ends.infercnv_obj"))
	keep <- unlist(lapply(unique(chr), function(ch) {
		i <- which(chr == ch)
		i[(5 + 1):(length(i) - 5)]
	}))
	expect_equal(rownames(step13@expr.data), genes[keep])
	expect_equal(step13@expr.data, step12@expr.data[genes[keep], ])
	## outlier pruning clips to the mean per-cell min and max of the object it is given
	step14 <- readRDS(file.path(od, "14_invert_log_transform.infercnv_obj"))
	step16 <- readRDS(file.path(od, "16_remove_outlier.infercnv_obj"))
	prelim <- readRDS(file.path(od, "preliminary.infercnv_obj"))
	expect_equal(prelim@expr.data, step14@expr.data)
	b <- infercnv:::.get_average_bounds(step14@expr.data)
	expect_equal(step16@expr.data, pmin(pmax(step14@expr.data, b[1]), b[2]))
	## masking replaces non-DE entries by a single value and leaves the rest alone
	step21 <- readRDS(file.path(od, "21_mask_nonDE.infercnv_obj"))
	changed <- step21@expr.data != step16@expr.data
	if (any(changed)) {
		expect_equal(length(unique(step21@expr.data[changed])), 1)
		expect_equal(unique(step21@expr.data[changed]), mean(step16@expr.data))
	}
	## reference cells are never masked
	ref_idx <- unlist(obj@reference_grouped_cell_indices)
	expect_equal(step21@expr.data[, ref_idx], step16@expr.data[, ref_idx])
	## denoising with a fixed filter sets values within 0.05 of the reference mean to that mean
	centre <- mean(step21@expr.data[, ref_idx])
	near <- abs(step21@expr.data - centre) < 0.05
	expect_true(all(res@expr.data[near] == centre))
	expect_equal(res@expr.data[!near], step21@expr.data[!near])
})
