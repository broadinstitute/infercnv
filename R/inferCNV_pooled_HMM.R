
# Pooled-count i6 HMM (experimental; run(HMM_expr_model = "pooled")).
#
# The default HMM works on the subcluster mean of per-cell log expression. Single-cell counts are mostly 0, 1 or 2,
# where log(x + 1) is compressed, so a one-copy gain moves a gene far less than log(1.5), and low-expressed genes much
# less than high-expressed ones. Here, per tumour subcluster:
#   1. counts (depth-normalised) are pooled over the subcluster's cells and compared with the pooled reference on the
#      log scale: a one-copy change moves every gene by ~log(1.5) or log(0.5), whatever its expression;
#   2. each gene's noise for a group of n cells is estimated from the reference cells: sampling (between-cell squared
#      coefficient of variation / n) plus, with several reference groups, the gene's variation between reference cell
#      types (which does not average out); both are moderated towards their trend on mean expression;
#   3. each gene's log ratio is limited to clip_z (3) sd, then the log ratios are smoothed along each chromosome with a small (pyramid) window, weighting genes by 1/sd
#      ('zscore': the mean of the genes' z-scores, rescaled to log-fold units), 1/variance ('invvar') or equally;
#   4. an i6 HMM is run with state means fixed at log(0.01), log(0.5), 0, log(1.5), log(2), log(3) and a
#      position-specific sd (the smoothed noise), scaled by a factor (log-linear in group size) calibrated on random
#      groups of reference cells (left out of the reference profile), so that their smoothed profiles are N(0, 1) in sd units; Student-t
#      emissions limit the pull of outlier positions.
# The same random reference groups give the null false-call rate.

.pooled_pyramid <- function(window_length) {
    h <- (window_length + 1) %/% 2
    c(seq_len(h), rev(seq_len(h - 1)))
}

# column-wise convolution with a symmetric kernel, zero-padded (weighted sums are then truncated at chromosome ends)
.pooled_conv <- function(M, k) {
    M <- as.matrix(M)
    h <- (length(k) - 1) %/% 2
    if (h == 0) return(M * k[1])
    P <- rbind(matrix(0, h, ncol(M)), M, matrix(0, h, ncol(M)))
    S <- stats::filter(P, k, sides = 2)
    S <- matrix(as.numeric(S), ncol = ncol(M))
    S[(h + 1):(h + nrow(M)), , drop = FALSE]
}

# Viterbi in log space: LL (positions x states) log-likelihoods
.pooled_viterbi <- function(LL, logT, logdelta) {
    n <- nrow(LL); K <- ncol(LL)
    V <- logdelta + LL[1, ]
    B <- matrix(0L, n, K)
    if (n > 1) {
        for (i in 2:n) {
            M <- V + logT                       # [from, to]
            b <- max.col(t(M), ties.method = "first")
            B[i, ] <- b
            V <- M[cbind(b, seq_len(K))] + LL[i, ]
        }
    }
    path <- integer(n)
    path[n] <- which.max(V)
    if (n > 1) for (i in (n - 1):1) path[i] <- B[i + 1, path[i + 1]]
    path
}

# moderate a positive per-gene quantity towards its loess trend on log mean expression (on the log scale)
.pooled_moderate <- function(x, mu, lambda) {
    ok <- is.finite(x) & x > 0 & mu > 0
    out <- x
    if (sum(ok) < 50) {
        out[!ok] <- stats::median(x[ok])
        return(out)
    }
    lx <- log(x[ok]); lm_ <- log(mu[ok])
    fit <- stats::loess(lx ~ lm_, span = 0.3, degree = 1)
    tr <- exp(stats::approx(lm_, stats::fitted(fit), xout = log(pmax(mu, min(mu[ok]))), rule = 2, ties = mean)$y)
    out[ok] <- exp((1 - lambda) * log(x[ok]) + lambda * log(tr[ok]))
    out[!ok] <- tr[!ok]
    out
}

# reference model: per-gene pooled reference profile and noise components
.pooled_reference_model <- function(N, ref_groups, moderation) {
    G <- length(ref_groups)
    m <- sapply(ref_groups, function(ix) Matrix::rowMeans(N[, ix, drop = FALSE]))
    m <- matrix(m, ncol = G)
    v <- sapply(ref_groups, function(ix) Matrix::rowMeans(N[, ix, drop = FALSE]^2))
    v <- matrix(v, ncol = G) - m^2
    ng <- lengths(ref_groups)
    mu <- rowMeans(m)
    v_within <- as.numeric(v %*% (ng / sum(ng)))
    cv2 <- ifelse(mu > 0, v_within / mu^2, NA)
    cv2 <- .pooled_moderate(cv2, mu, moderation)
    b <- rep(0, length(mu))
    if (G >= 2) {
        lm_ <- log(sweep(m, 2, 0.5 / ng, "+"))
        b <- apply(lm_, 1, stats::var)
        b <- .pooled_moderate(b, mu, moderation)
    }
    list(m = m, ng = ng, mu = mu, cv2 = cv2, b = b, n_ref = sum(ng))
}

# smoothed log-ratio profiles and their sd for groups of cells (columns: groups); mu_list optionally gives a per-group
# reference profile (leave-out)
.pooled_profiles <- function(N, groups, model, chr_idx, window_length, gene_weights, center, mu_list = NULL, clip_z = 3) {
    k <- .pooled_pyramid(window_length)
    S <- matrix(NA_real_, nrow(N), length(groups)); SD <- S
    for (j in seq_along(groups)) {
        ix <- groups[[j]]; n <- length(ix)
        mu <- if (is.null(mu_list)) model$mu else mu_list[[j]]
        msub <- Matrix::rowMeans(N[, ix, drop = FALSE])
        r <- log((msub + 0.5 / n) / (mu + 0.5 / model$n_ref))
        s2 <- model$cv2 / n + model$b
        w <- switch(gene_weights, zscore = 1 / sqrt(s2), invvar = 1 / s2, equal = rep(1, length(s2)))
        if (center) {
            reliable <- mu * n >= 5
            if (sum(reliable) < 100) reliable <- rep(TRUE, length(r))
            r <- r - stats::median(r[reliable])
        }
        ## limit single genes (co-expressed gene clusters, outlier cells) to clip_z sd; a one-copy change is < 1 sd per gene
        if (is.finite(clip_z)) {
            lim <- clip_z * sqrt(s2)
            r <- pmin(pmax(r, -lim), lim)
        }
        for (gi in chr_idx) {
            num <- .pooled_conv(w[gi] * r[gi], k)
            den <- .pooled_conv(w[gi], k)
            vv <- .pooled_conv(w[gi]^2 * s2[gi], k^2)
            S[gi, j] <- num / den
            SD[gi, j] <- sqrt(vv) / den
        }
    }
    list(S = S, SD = SD)
}

.pooled_states <- function(S, SD, kappa, chr_idx, t, df) {
    lev <- log(c(0.01, 0.5, 1, 1.5, 2, 3))
    logT <- log(matrix(t, 6, 6) + diag(1 - 6 * t, 6))
    logdelta <- log(c(t, t, 1 - 5 * t, t, t, t))
    st <- matrix(NA_integer_, nrow(S), ncol(S))
    for (j in seq_len(ncol(S))) {
        for (gi in chr_idx) {
            sc <- kappa[j] * SD[gi, j]
            LL <- sapply(lev, function(m) stats::dt((S[gi, j] - m) / sc, df = df, log = TRUE) - log(sc))
            st[gi, j] <- .pooled_viterbi(matrix(LL, ncol = 6), logT, logdelta)
        }
    }
    st
}

predict_CNV_via_pooled_HMM_on_tumor_subclusters <- function(infercnv_obj, window_length = 31,
                                                            gene_weights = c("zscore", "invvar", "equal"),
                                                            t = 1e-6, moderation = 0.5, df = 5, center = TRUE,
                                                            max_null_groups = 60, min_null_groups = 5) {
    gene_weights <- match.arg(gene_weights)
    flog.info(sprintf("pooled-count i6 HMM: window %d genes, %s gene weights", window_length, gene_weights))
    if (is.null(infercnv_obj@tumor_subclusters)) stop("pooled HMM needs tumour subclusters (analysis_mode='subclusters')")
    C <- infercnv_obj@count.data
    if (!identical(rownames(C), rownames(infercnv_obj@expr.data))) C <- C[rownames(infercnv_obj@expr.data), colnames(infercnv_obj@expr.data)]
    cs <- Matrix::colSums(C)
    N <- C %*% Matrix::Diagonal(x = stats::median(cs) / pmax(cs, 1))
    colnames(N) <- colnames(C)
    ref_groups <- infercnv_obj@reference_grouped_cell_indices
    if (length(ref_groups) == 0) stop("pooled HMM needs reference cells")
    model <- .pooled_reference_model(N, ref_groups, moderation)
    gene_order <- infercnv_obj@gene_order
    chr_idx <- split(seq_len(nrow(gene_order)), factor(gene_order$chr, levels = unique(gene_order$chr)))

    ## null groups of reference cells, sized like the observation subclusters, left out of the reference profile
    sub <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive = FALSE)
    ref_cells <- unlist(ref_groups)
    sizes <- lengths(sub[!vapply(sub, function(ix) all(ix %in% ref_cells), logical(1))])
    if (exists(".Random.seed", envir = globalenv())) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    }
    set.seed(1)
    draw <- function(sizes) {
        groups <- list(); src <- integer(0)
        for (g in seq_along(ref_groups)) {
            pool <- ref_groups[[g]]; pool <- pool[sample.int(length(pool))]
            ## leave at least half of each reference group in the reference profile
            avail <- length(pool) %/% 2
            used <- 0
            while (avail - used >= min(sizes) && length(groups) < max_null_groups) {
                k <- min(sizes[sample.int(length(sizes), 1)], avail - used)
                groups[[length(groups) + 1]] <- pool[used + seq_len(k)]; src <- c(src, g)
                used <- used + k
            }
        }
        list(groups = groups, src = src)
    }
    nd <- draw(sizes)
    if (length(nd$groups) < min_null_groups) {
        cap <- max(floor(sum(lengths(ref_groups)) / (4 * min_null_groups)), 2)
        nd <- draw(pmin(sizes, cap))
    }
    if (length(nd$groups) < min_null_groups) stop("pooled HMM: too few reference cells for the null calibration")
    G <- length(ref_groups)
    mu_list <- lapply(seq_along(nd$groups), function(j) {
        g <- nd$src[j]; ix <- nd$groups[[j]]
        mg <- (model$ng[g] * model$m[, g] - length(ix) * Matrix::rowMeans(N[, ix, drop = FALSE])) / (model$ng[g] - length(ix))
        model$mu + (mg - model$m[, g]) / G
    })
    np <- .pooled_profiles(N, nd$groups, model, chr_idx, window_length, gene_weights, center, mu_list)
    ## noise scale per group size: spread of each null profile in model sd units, log-linear in log group size
    ## (correlated neighbouring genes inflate the smoothed noise by a size-dependent factor)
    kj <- vapply(seq_along(nd$groups), function(j) { z <- np$S[, j] / np$SD[, j]; stats::mad(z[is.finite(z)], center = 0) },
                 numeric(1))
    nj <- lengths(nd$groups)
    kappa_at <- if (length(unique(nj)) >= 3 && max(nj) / min(nj) >= 1.5) {
        fit <- stats::lm(log(kj) ~ log(nj))
        function(n) pmin(pmax(exp(stats::predict(fit, newdata = data.frame(nj = n))), min(kj)), max(kj))
    } else {
        function(n) rep(stats::median(kj), length(n))
    }
    null_st <- .pooled_states(np$S, np$SD, kappa_at(nj), chr_idx, t, df)
    null_fpr <- mean(null_st != 3)
    flog.info(sprintf(paste0("pooled HMM: noise scale %.2f-%.2f (groups of %d-%d cells; from %d reference groups); ",
                             "null false-call rate %.1f%%"), min(kappa_at(range(nj))), max(kappa_at(range(nj))), min(nj), max(nj),
                      length(nd$groups), 100 * null_fpr))

    ## all subclusters (observation and reference)
    op <- .pooled_profiles(N, sub, model, chr_idx, window_length, gene_weights, center)
    kappa <- kappa_at(lengths(sub))
    st <- .pooled_states(op$S, op$SD, kappa, chr_idx, t, df)
    hmm.data <- infercnv_obj@expr.data
    hmm.data[, ] <- -1
    for (j in seq_along(sub)) hmm.data[, sub[[j]]] <- st[, j]
    infercnv_obj@expr.data <- hmm.data
    infercnv_obj@options$pooled_hmm <- list(window_length = window_length, gene_weights = gene_weights, kappa = kappa,
                                            null_fpr = null_fpr, null_groups = length(nd$groups),
                                            profiles = op$S, profile_sd = op$SD)
    infercnv_obj
}
