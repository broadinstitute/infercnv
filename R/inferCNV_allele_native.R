# Bayesian models for the allele-aware steps (step 24, allele-only filter; step 27, combined expression + allele model),
# sampled in R.
#
# Each CNV region (a set of genes in one group of cells, as reported by the HMMs) is modelled as a mixture over
# copy-number states. Cell j has a state eps_j ~ Cat(theta), theta ~ Dirichlet(1, ..., 1), and each state is a
# copy-number genotype with an expected lesser-allele fraction f:
#
#   i6 (+ cnLOH)   1: 0 copies (allele data uninformative, treated as balanced)   2: 1 copy (1,0), f = err
#                  3: 2 copies (1,1), f = 1/2      4: 3 copies (2,1), f = 1/3     5: 4 copies (2,2) or (3,1), f = 1/2 or 1/4
#                  6: >= 5 copies, f = 1/2, 1/3 or 1/6      7: copy-neutral LOH (2,0), f = err
#   i3 (+ cnLOH)   1: loss, f = err   2: neutral, f = 1/2   3: gain, f = 1/2, 1/3 or 1/4   4: copy-neutral LOH, f = err
#   allele only    1: allelic imbalance, f = err, 1/3 or 1/4   2: balanced, f = 1/2
#
# (states with several genotypes give each the same weight). Allele counts: for SNP s in cell j, the lesser-allele count
# r_sj of n_sj reads is beta-binomial with mean p and concentration kappa, where p = f if the lesser allele of SNP s is
# the one present in fewer copies in this region (a_s = 1) and p = 1 - f otherwise; a_s ~ Bernoulli(1/2) is shared by
# all cells of the region, so a loss of heterozygosity must keep the same allele across cells. kappa (allelic
# overdispersion, including transcriptional bursting) is fitted by maximum likelihood on the reference cells, where
# SNPs are balanced. Expression: the cell's mean residual over the region's genes, x_j, is normal with the state's
# h-spike mean mu_k and standard deviation s_ref * sigma_k / sigma_neutral, where s_ref is the standard deviation of
# the same region mean across reference cells. Using the region mean, with its spread measured on reference cells,
# accounts for the correlation of neighbouring genes after smoothing, so that expression and allele evidence are
# weighted by how informative each actually is.
#
# Gibbs sampling alternates (eps_j, genotype component of cell j) | theta, a, data; theta | eps; and a_s | eps, data.
# All per-SNP terms are precomputed once, so a sweep costs O(SNP-cell pairs with coverage). Output is an mcmc.list
# with columns theta[k] and epsilon[j], as rjags::coda.samples() returns, so that the downstream code is unchanged.

.allele_state_specs <- function(type = c("i6", "i3", "i2"), enable_cnLOH = TRUE, err = 0.01) {
    type <- match.arg(type)
    s <- switch(type,
        i6 = list(list(f = 0.5), list(f = err), list(f = 0.5), list(f = 1/3), list(f = c(0.5, 0.25)),
                  list(f = c(0.5, 1/3, 1/6)), list(f = err)),
        i3 = list(list(f = err), list(f = 0.5), list(f = c(0.5, 1/3, 0.25)), list(f = err)),
        i2 = list(list(f = c(err, 1/3, 0.25)), list(f = 0.5)))
    if (type == "i6" && !enable_cnLOH) s <- s[1:6]
    if (type == "i3" && !enable_cnLOH) s <- s[1:3]
    s
}

# beta-binomial log-likelihood without the binomial coefficient (the same for every state)
.bb_ll <- function(r, n, p, kappa) {
    lbeta(r + kappa * p, n - r + kappa * (1 - p)) - lbeta(kappa * p, kappa * (1 - p))
}

# Maximum-likelihood concentration of a beta-binomial over the SNP-cell pairs with coverage: with mean 1/2 (reference
# cells, where heterozygous SNPs are balanced), or with each SNP's mean set to its pooled lesser-allele fraction
# (free_mean = TRUE; used when no reference cell has allele data, so that allelic imbalance from CNAs is not
# mistaken for overdispersion).
.fit_allele_kappa <- function(R, N, free_mean = FALSE) {
    N <- as(N, "CsparseMatrix"); tr <- Matrix::summary(N)
    if (nrow(tr) < 50) return(10)
    r <- as.numeric(R[cbind(tr$i, tr$j)]); n <- tr$x
    p <- if (free_mean) {
        ps <- rowsum(r, tr$i)[, 1] / rowsum(n, tr$i)[, 1]
        pmin(pmax(ps[as.character(tr$i)], 0.01), 0.99)
    } else 0.5
    f <- function(lk) -sum(.bb_ll(r, n, p, exp(lk)))
    exp(stats::optimize(f, c(log(0.05), log(1000)))$minimum)
}

# Per-region allele terms: for f = 1/2 a per-cell sum; for every other f the per-cell sum with a_s = 0 (p = 1 - f)
# and the sparse SNP x cell difference (a_s = 1 minus a_s = 0).
.allele_region_terms <- function(R, N, fracs, kappa) {
    N <- as(N, "CsparseMatrix"); tr <- Matrix::summary(N)
    S <- nrow(N); C <- ncol(N)
    r <- as.numeric(R[cbind(tr$i, tr$j)]); n <- tr$x
    out <- list()
    for (f in unique(fracs)) {
        key <- sprintf("%.10g", f)
        if (isTRUE(all.equal(f, 0.5))) {
            out[[key]] <- list(base = .tabulate_sum(tr$j, .bb_ll(r, n, 0.5, kappa), C), D = NULL)
        } else {
            l1 <- .bb_ll(r, n, f, kappa); l0 <- .bb_ll(r, n, 1 - f, kappa)
            out[[key]] <- list(base = .tabulate_sum(tr$j, l0, C),
                               D = Matrix::sparseMatrix(i = tr$i, j = tr$j, x = l1 - l0, dims = c(S, C)))
        }
    }
    out
}

.tabulate_sum <- function(idx, x, n) {
    v <- numeric(n)
    s <- rowsum(x, idx)
    v[as.integer(rownames(s))] <- s[, 1]
    v
}

# Gibbs sampler. expr_ll: K x C expression log-likelihoods (zeros for the allele-only model); allele: output of
# .allele_region_terms or NULL; states: .allele_state_specs(). One chain per state, initialised with all cells in
# that state and every a_s = 1.
.native_allele_gibbs <- function(expr_ll, allele, states, n_burn = 700, n_keep = 1000) {
    K <- length(states); C <- ncol(expr_ll); M <- K
    # flatten (state, genotype component) pairs
    q_state <- unlist(lapply(seq_len(K), function(k) rep(k, length(states[[k]]$f))))
    q_f <- unlist(lapply(states, `[[`, "f"))
    q_lw <- unlist(lapply(states, function(s) rep(-log(length(s$f)), length(s$f))))
    q_key <- sprintf("%.10g", q_f)
    Q <- length(q_state)
    has_allele <- !is.null(allele)
    imb_keys <- if (has_allele) names(allele)[!vapply(allele, function(a) is.null(a$D), logical(1))] else character(0)
    S <- if (length(imb_keys)) nrow(allele[[imb_keys[1]]]$D) else 0L
    chain <- rep(seq_len(M), each = C)
    a <- matrix(1, max(S, 1L), M)                    # SNP x chain
    eps <- chain                                       # chain m starts with all cells in state m
    q_cur <- match(eps, q_state)                       # first component of that state
    base_q <- expr_ll[q_state, , drop = FALSE] + q_lw  # Q x C
    theta_draws <- array(0, c(n_keep, K, M)); eps_draws <- array(0L, c(n_keep, C, M))
    for (it in seq_len(n_burn + n_keep)) {
        cnt <- matrix(tabulate(eps + K * (chain - 1L), K * M), K, M)
        g <- matrix(stats::rgamma(K * M, 1 + cnt), K, M); th <- g / rep(colSums(g), each = K)
        # log-likelihood of each (state, component) for every cell and chain: Q x (C * M)
        LL <- base_q[, rep(seq_len(C), M), drop = FALSE]
        if (has_allele) {
            for (key in unique(q_key)) {
                al <- allele[[key]]
                v <- if (is.null(al$D)) rep(al$base, M) else
                     as.vector(al$base + as.matrix(Matrix::crossprod(al$D, a)))
                rows <- which(q_key == key)
                LL[rows, ] <- LL[rows, , drop = FALSE] + rep(v, each = length(rows))
            }
        }
        LL <- LL + log(th[q_state, chain, drop = FALSE])
        LL <- LL - rep(apply(LL, 2, max), each = Q)
        P <- exp(LL)
        cp <- P; for (q in seq_len(Q)[-1]) cp[q, ] <- cp[q - 1, ] + P[q, ]
        u <- stats::runif(C * M) * cp[Q, ]
        q_cur <- 1L + colSums(cp < rep(u, each = Q))
        eps <- q_state[q_cur]
        # direction of imbalance per SNP, given the cells' genotype components
        if (S > 0) {
            logit <- matrix(0, S, M)
            for (key in imb_keys) {
                ind <- matrix(as.numeric(q_key[q_cur] == key), C, M)
                logit <- logit + as.matrix(allele[[key]]$D %*% ind)
            }
            a <- matrix(as.numeric(stats::runif(S * M) < stats::plogis(logit)), S, M)
        }
        if (it > n_burn) {
            theta_draws[it - n_burn, , ] <- th
            eps_draws[it - n_burn, , ] <- matrix(eps, C, M)
        }
    }
    cols <- c(sprintf("theta[%d]", seq_len(K)), sprintf("epsilon[%d]", seq_len(C)))
    coda::mcmc.list(lapply(seq_len(M), function(m) {
        coda::mcmc(matrix(c(theta_draws[, , m], eps_draws[, , m]), n_keep, K + C, dimnames = list(NULL, cols)),
                   start = n_burn + 1)
    }))
}

# Expression log-likelihoods (K x C) of the per-cell region means, with the spread of the region mean measured on the
# reference cells and scaled per state by the h-spike single-cell sd ratio. sig: precisions (1/sd^2) of the states;
# extra states beyond length(mu) (cnLOH) use the neutral mean and sd.
.region_expr_ll <- function(expr_obs, expr_ref, mu, sig, neutral, K) {
    if (nrow(expr_obs) == 0) return(matrix(0, K, ncol(expr_obs)))   # no genes: allele evidence only
    xbar <- colMeans(as.matrix(expr_obs))
    sd1 <- 1 / sqrt(sig)
    s_ref <- if (!is.null(expr_ref) && ncol(expr_ref) >= 5) stats::sd(colMeans(as.matrix(expr_ref))) else sd1[neutral]
    if (!is.finite(s_ref) || s_ref <= 0) s_ref <- sd1[neutral]
    mu_k <- c(mu, rep(mu[neutral], K - length(mu)))
    sd_k <- s_ref * c(sd1, rep(sd1[neutral], K - length(mu))) / sd1[neutral]
    # K x C, also when the region has a single cell
    do.call(rbind, lapply(seq_len(K), function(k) stats::dnorm(xbar, mu_k[k], sd_k[k], log = TRUE)))
}

# reference cells by name, and kappa fitted on their allele counts
.allele_reference <- function(infercnv_allele_obj) {
    ref_cells <- colnames(infercnv_allele_obj@expr.data)[unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
    a_ref <- intersect(ref_cells, colnames(infercnv_allele_obj@allele.count.data))
    free_mean <- length(a_ref) == 0 ||
        Matrix::nnzero(infercnv_allele_obj@allele.coverage.data[, a_ref, drop = FALSE]) < 50
    if (free_mean) {
        flog.warn("Reference cells have no allele data; fitting allelic overdispersion on all cells with per-SNP means")
        a_ref <- colnames(infercnv_allele_obj@allele.count.data)
    }
    kappa <- .fit_allele_kappa(infercnv_allele_obj@allele.count.data[, a_ref, drop = FALSE],
                               infercnv_allele_obj@allele.coverage.data[, a_ref, drop = FALSE], free_mean)
    flog.info(sprintf("Allelic overdispersion fitted on %d cells: beta-binomial concentration %.3g", length(a_ref), kappa))
    list(cells = ref_cells, kappa = kappa)
}

# mclapply returns try-error objects for failed regions; stop with the first message instead of failing downstream
.check_mc <- function(res) {
    bad <- vapply(res, function(r) inherits(r, "try-error") || is.null(r), logical(1))
    if (any(bad)) {
        msg <- if (inherits(res[[which(bad)[1]]], "try-error")) as.character(res[[which(bad)[1]]]) else "no result returned"
        stop(sprintf("%d of %d regions failed in the native sampler; first error: %s", sum(bad), length(res), msg))
    }
    res
}

# Step 24: allele-only model per region (states: 1 imbalanced, 2 balanced).
.run_allele_native <- function(MCMC_inferCNV_obj, infercnv_allele_obj, allele_mode, cores) {
    ref <- .allele_reference(infercnv_allele_obj)
    states <- .allele_state_specs("i2")
    fr <- unlist(lapply(states, `[[`, "f"))
    .check_mc(parallel::mclapply(MCMC_inferCNV_obj@cell_gene, function(x) {
        rows <- if (allele_mode == "snp_level") x$SNPs else x$Genes
        cells <- x$Cells
        al <- .allele_region_terms(infercnv_allele_obj@allele.count.data[rows, cells, drop = FALSE],
                                   infercnv_allele_obj@allele.coverage.data[rows, cells, drop = FALSE], fr, ref$kappa)
        .native_allele_gibbs(matrix(0, length(states), length(cells)), al, states)
    }, mc.cores = cores))
}

# Step 27: combined expression + allele model per region.
.run_combined_native <- function(MCMC_inferCNV_combined_obj, infercnv_allele_obj, HMM_type, enable_cnLOH, cores) {
    ref <- .allele_reference(infercnv_allele_obj)
    expr <- infercnv_allele_obj@expr.data
    ref_idx <- match(ref$cells, colnames(expr))
    mu <- MCMC_inferCNV_combined_obj@infercnv_obj_mu
    sig <- MCMC_inferCNV_combined_obj@infercnv_obj_sig
    neutral <- if (HMM_type == "i6") 3 else 2
    .check_mc(parallel::mclapply(MCMC_inferCNV_combined_obj@cell_gene, function(x) {
        rows <- if (!is.null(x$infercnv_allele_SNPs)) x$infercnv_allele_SNPs else x$infercnv_allele_Genes
        has_allele <- length(rows) > 0
        states <- .allele_state_specs(HMM_type, enable_cnLOH && has_allele)
        cells <- colnames(expr)[x$infercnv_Cells]
        expr_ll <- .region_expr_ll(expr[x$infercnv_Genes, x$infercnv_Cells, drop = FALSE],
                                   if (length(ref_idx)) expr[x$infercnv_Genes, ref_idx, drop = FALSE] else NULL,
                                   mu, sig, neutral, length(states))
        al <- NULL
        if (has_allele) {
            acells <- match(cells, colnames(infercnv_allele_obj@allele.count.data))
            stopifnot(!anyNA(acells))
            al <- .allele_region_terms(infercnv_allele_obj@allele.count.data[rows, acells, drop = FALSE],
                                       infercnv_allele_obj@allele.coverage.data[rows, acells, drop = FALSE],
                                       unlist(lapply(states, `[[`, "f")), ref$kappa)
        }
        .native_allele_gibbs(expr_ll, al, states)
    }, mc.cores = cores))
}
