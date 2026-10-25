# Clones from subcluster HMM calls, with CNV events that share their boundaries across clones.
#
# inferCNV calls CNVs per tumour subcluster, and each subcluster's HMM places its own boundaries, so one event
# inherited by several subclusters is reported with slightly different edges in each. Starting from those calls,
#
#  1. Shared boundaries. Changes of direction (loss / neutral / gain) in the subclusters' HMM calls that fall within
#     `tol` genes of each other are taken to be one breakpoint, placed once where it best separates the mean
#     unsmoothed residuals of all the subclusters that carry it; each subcluster is then re-called on the shared
#     segments (most common direction, then most common state).
#  2. Candidate clones. Subclusters with the same calls form a candidate clone; two candidates are merged unless some
#     segment differs in direction with a clear difference in the cells' residuals (Welch t and a fraction of the
#     one-copy step).
#  3. Parsimony re-call. Each clone is re-called on the shared segments by a Viterbi pass along each chromosome that
#     minimises the squared distance between the clone's segment level and the level of the chosen state (weighted by
#     the number of genes) plus `lambda` per change of state; allowed states at a segment are neutral and the states
#     any clone was called with there. A per-segment offset common to all clones (expression of the region's genes in
#     these cells, not copy number) is removed first. Clones whose re-calls differ only by events of fewer than
#     `min_event_genes` genes, or by less than `d_frac` of a one-copy step (median per-gene difference), are merged
#     and the re-call repeated.
#  4. Clone tree. Each event (consecutive segments with the same direction and the same carrier clones) has a carrier
#     set; nested carrier sets give the tree (events carried by all clones are truncal, the others derived). Carrier
#     sets that conflict with the tree are reported as unplaced.
#
# States are coded as i6 (1-6, neutral 3); i3 calls (1-3) are mapped to 2, 3, 4 and back on output.


#' @title infercnv_clones
#'
#' @description Group tumour subclusters into clones by shared CNV events, refine event boundaries so that an event
#' has the same boundaries in every clone that carries it, re-call each clone with a parsimony penalty, and build the
#' clone tree (truncal and derived events). Runs on the HMM calls of a subcluster-mode run.
#'
#' @param infercnv_obj infercnv object with subclusters, raw counts (count.data) and reference groups (e.g. the
#'                     object at the HMM step of run(), or run.final.infercnv_obj)
#' @param hmm_states genes x cells matrix of HMM states (i6: 1-6, i3: 1-3), e.g. the expr.data of the HMM object
#'                   (after the Bayesian filter when used)
#' @param HMM_type "i6" or "i3"
#' @param out_dir directory for the output files (NULL: no files)
#' @param output_prefix prefix of the output files
#' @param pool_groups optional list of character vectors: observation groups analysed together (e.g. the epithelial
#'                    clusters of a spatial sample); other groups are analysed separately
#' @param min_cells subclusters with fewer cells do not form clones of their own and join the nearest clone (default: 5)
#' @param tol breakpoints of different subclusters within this many genes are candidates for one shared breakpoint (default: 30)
#' @param lambda penalty per change of state in the parsimony re-call (default: 10)
#' @param min_event_genes smallest event (genes) that can separate two clones (default: 100)
#' @param d_frac smallest difference, as a fraction of the one-copy step, that can separate two clones (default: 0.7)
#' @param plot write a clone tree + clone copy-number figure per group (default: TRUE)
#'
#' @return list per group: cells (cell, group, subcluster, clone), states (genes x clones, after the parsimony
#'         re-call, in the input HMM coding), events, tree, unplaced events
#'
#' @export
#'
infercnv_clones <- function(infercnv_obj, hmm_states, HMM_type = c("i6", "i3"), out_dir = NULL,
                            output_prefix = "infercnv.clones", pool_groups = NULL, min_cells = 5, tol = 30,
                            lambda = 10, min_event_genes = 100, d_frac = 0.7, plot = TRUE) {
    HMM_type <- match.arg(HMM_type)
    genes <- intersect(rownames(hmm_states), rownames(infercnv_obj@expr.data))
    gene_order <- infercnv_obj@gene_order[genes, , drop = FALSE]
    chr <- as.character(gene_order$chr); pos <- gene_order$start
    cells <- colnames(infercnv_obj@expr.data)
    hmm <- as.matrix(hmm_states[genes, cells, drop = FALSE])
    if (HMM_type == "i3") hmm <- hmm + 1
    R <- .clone_unsmoothed_residuals(infercnv_obj, genes)

    sc <- .clone_subclusters(infercnv_obj)
    if (!is.null(pool_groups)) for (k in seq_along(pool_groups)) {
        nm <- paste(pool_groups[[k]], collapse = "+")
        sel <- sc$group %in% pool_groups[[k]]
        sc$sub[sel] <- paste(sc$group[sel], sc$sub[sel], sep = "|"); sc$group[sel] <- nm
    }
    res <- list()
    for (g in unique(sc$group)) {
        s_g <- sc[sc$group == g, ]
        subs <- unique(s_g$sub); cbs <- split(s_g$cell, factor(s_g$sub, levels = subs)); ns <- lengths(cbs)
        P <- vapply(subs, function(s) hmm[, cbs[[s]][1]], numeric(length(genes)))
        if (!is.matrix(P)) P <- matrix(P, ncol = 1, dimnames = list(NULL, subs))
        M <- vapply(subs, function(s) rowMeans(R[, cbs[[s]], drop = FALSE]), numeric(length(genes)))
        if (!is.matrix(M)) M <- matrix(M, ncol = 1, dimnames = list(NULL, subs))
        keep <- subs[ns >= 3]
        if (!length(keep)) next
        flog.info(sprintf("infercnv_clones: %s, %d cells in %d subclusters", g, sum(ns), length(subs)))
        rf <- .clone_refine_group(P[, keep, drop = FALSE], M[, keep, drop = FALSE], chr, tol = tol)
        P2 <- P; P2[, keep] <- rf$P
        seg <- rf$seg
        cl <- if (sum(ns >= min_cells) >= 2) .clone_candidates(P2, seg, R, cbs, min_cells = min_cells)
              else list(clone = setNames(rep("A", length(subs)), subs), low = character(0))
        mp <- .clone_merge(cl, P2, seg, chr, R, cbs, lambda = lambda, min_event_genes = min_event_genes, d_frac = d_frac)
        ev <- .clone_shared_events(mp$pcons, seg, chr)
        ev$start_pos <- pos[ev$start]; ev$end_pos <- gene_order$stop[ev$end]
        ev$start_gene <- genes[ev$start]; ev$end_gene <- genes[ev$end]
        tl <- .clone_tree_layout(ev, mp$sizes)
        states <- mp$pcons; if (HMM_type == "i3") states <- states - 1
        rownames(states) <- genes
        cells_df <- data.frame(cell = cells[s_g$cell], group = g, subcluster = s_g$sub,
                               clone = mp$clone[s_g$sub], stringsAsFactors = FALSE)
        flog.info(sprintf("infercnv_clones: %s -> %d clones (%s); %d events, %d not placed on the tree", g,
                          length(mp$sizes), paste(names(mp$sizes), mp$sizes, sep = "=", collapse = " "),
                          nrow(ev), nrow(tl$unplaced)))
        res[[g]] <- list(cells = cells_df, states = states, events = ev, tree = tl$tree, unplaced = tl$unplaced,
                         breakpoints = rf$breaks, sizes = mp$sizes)
    }

    if (!is.null(out_dir) && length(res)) {
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        f <- function(x) file.path(out_dir, paste0(output_prefix, x))
        cells_out <- do.call(rbind, lapply(res, `[[`, "cells"))
        write.table(cells_out, f(".cells.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
        ev_out <- do.call(rbind, lapply(names(res), function(g) { e <- res[[g]]$events
            if (!nrow(e)) return(NULL)
            data.frame(group = g, chr = e$chr, start = e$start_pos, end = e$end_pos, start_gene = e$start_gene,
                       end_gene = e$end_gene, genes = e$genes, direction = ifelse(e$dir > 0, "gain", "loss"),
                       clones = e$carriers, truncal = e$carriers == paste(sort(names(res[[g]]$sizes)), collapse = ""),
                       stringsAsFactors = FALSE) }))
        if (!is.null(ev_out)) write.table(ev_out, f(".events.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
        tr_out <- do.call(rbind, lapply(names(res), function(g) { t <- res[[g]]$tree
            data.frame(group = g, node = t$node, parent = t$parent, cells = t$size, genes_changed = t$len,
                       events = t$events, stringsAsFactors = FALSE) }))
        write.table(tr_out, f(".tree.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
        st <- do.call(cbind, lapply(names(res), function(g) { s <- res[[g]]$states
            colnames(s) <- paste(g, colnames(s), sep = "."); s }))
        write.table(st, f(".states.dat"), sep = "\t", quote = FALSE)
        if (plot) for (g in names(res)) {
            if (length(res[[g]]$sizes) < 1) next
            fn <- f(sprintf(".%s.png", gsub("[^A-Za-z0-9_.+-]", "_", g)))
            .clone_plot(res[[g]], chr, HMM_type, title = g, file = fn)
        }
    }
    invisible(res)
}


## ---- inputs ----

# unsmoothed log residuals vs the reference: depth-normalised log2(x + 1) minus the mean of the reference group
# means, capped at +-cap and centred per cell (as inferCNV centres cells after smoothing)
.clone_unsmoothed_residuals <- function(infercnv_obj, genes, cap = 1.5) {
    cells <- colnames(infercnv_obj@expr.data)
    cs <- Matrix::colSums(infercnv_obj@count.data[, cells, drop = FALSE])
    cnt <- infercnv_obj@count.data[genes, cells, drop = FALSE]
    x <- as.matrix(log2(sweep(as.matrix(cnt), 2, cs, "/") * stats::median(cs) + 1))
    refm <- rowMeans(vapply(infercnv_obj@reference_grouped_cell_indices,
                            function(i) rowMeans(x[, i, drop = FALSE]), numeric(nrow(x))))
    r <- x - refm
    r[r > cap] <- cap; r[r < -cap] <- -cap
    sweep(r, 2, apply(r, 2, stats::median))
}

.clone_subclusters <- function(infercnv_obj) {
    do.call(rbind, lapply(names(infercnv_obj@observation_grouped_cell_indices), function(g) {
        sc <- infercnv_obj@tumor_subclusters$subclusters[[g]]
        if (is.null(sc)) sc <- list(all = infercnv_obj@observation_grouped_cell_indices[[g]])
        do.call(rbind, lapply(names(sc), function(s) data.frame(cell = as.integer(sc[[s]]), group = g, sub = s,
                                                                stringsAsFactors = FALSE)))
    }))
}

.clone_dir <- function(P) sign(P - 3)

.clone_mode <- function(x) { t <- table(x); as.numeric(names(t)[which.max(t)]) }

.clone_wmode <- function(x, w) { t <- tapply(w, x, sum); as.numeric(names(t)[which.max(t)]) }

.clone_wmedian <- function(x, w) { o <- order(x); x[o][which(cumsum(w[o]) >= sum(w) / 2)[1]] }


## ---- 1. shared breakpoints ----

# one row per change of state along a chromosome (pos = first gene of the new segment)
.clone_breaks <- function(P, chr) {
    out <- lapply(colnames(P), function(s) {
        v <- P[, s]
        i <- which(v[-1] != v[-length(v)] & chr[-1] == chr[-length(chr)]) + 1L
        if (!length(i)) return(NULL)
        data.frame(sub = s, chr = chr[i], pos = i, stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

# 1-D clustering of breakpoints per chromosome: gaps <= tol join; clusters wider than 3 * tol split at the largest gap
.clone_cluster_breaks <- function(bk, tol) {
    bk <- bk[order(bk$pos), ]
    split_rec <- function(p) {
        if (length(p) < 2 || diff(range(p)) <= 3 * tol) return(list(seq_along(p)))
        g <- which.max(diff(p))
        c(split_rec(p[1:g]), lapply(split_rec(p[(g + 1):length(p)]), function(i) i + g))
    }
    bk$cl <- NA_character_
    for (ch in unique(bk$chr)) {
        ii <- which(bk$chr == ch); p <- bk$pos[ii]
        g <- cumsum(c(1L, diff(p) > tol)); k <- 0L
        for (u in unique(g)) { jj <- which(g == u)
            for (pp in split_rec(p[jj])) { k <- k + 1L; bk$cl[ii[jj[pp]]] <- paste(ch, k, sep = ":") } }
    }
    bk
}

# change-point statistic at each candidate p for the columns of Y, using rows lo..hi-1
.clone_cp_stats <- function(Y, lo, hi, cand) {
    Z <- Y[lo:(hi - 1), , drop = FALSE]
    C <- apply(Z, 2, cumsum); if (!is.matrix(C)) C <- matrix(C, nrow = 1)
    tot <- C[nrow(C), ]
    S <- vapply(cand, function(p) { nl <- p - lo; nr <- hi - p
        ml <- C[nl, ] / nl; mr <- (tot - C[nl, ]) / nr; nl * nr / (nl + nr) * (ml - mr)^2 }, numeric(ncol(Y)))
    if (!is.matrix(S)) S <- matrix(S, nrow = ncol(Y))
    S
}

# place each breakpoint cluster once, on the residuals of all its carrier subclusters (weighted by their noise)
.clone_place_breaks <- function(bk, M, chr, flank = 150, min_seg = 5, tol_window = 30L) {
    sd_s <- apply(M, 2, function(y) stats::mad(diff(y), na.rm = TRUE) / sqrt(2)); w <- 1 / pmax(sd_s, 1e-6)^2
    cls <- do.call(rbind, lapply(split(bk, bk$cl), function(b)
        data.frame(cl = b$cl[1], chr = b$chr[1], init = as.integer(round(stats::median(b$pos))), wlo = min(b$pos),
                   whi = max(b$pos), stringsAsFactors = FALSE)))
    carriers <- lapply(split(bk$sub, bk$cl), unique)
    cls <- cls[order(cls$init), ]
    chr_lo <- tapply(seq_along(chr), chr, min); chr_hi <- tapply(seq_along(chr), chr, max) + 1L
    cls$pos <- cls$init; cls$n_carriers <- 0L; cls$spread <- NA_real_
    for (k in seq_len(nrow(cls))) {
        ch <- cls$chr[k]
        prv <- if (k > 1 && cls$chr[k - 1] == ch) cls$init[k - 1] else chr_lo[[ch]]
        nxt <- if (k < nrow(cls) && cls$chr[k + 1] == ch) cls$init[k + 1] else chr_hi[[ch]]
        lo <- max(prv, cls$wlo[k] - flank); hi <- min(nxt, cls$whi[k] + flank)
        car <- carriers[[cls$cl[k]]]; cls$n_carriers[k] <- length(car)
        a <- max(lo + min_seg, cls$wlo[k] - tol_window); b <- min(hi - min_seg, cls$whi[k] + tol_window)
        if (b < a || hi - lo < 2 * min_seg) next
        cand <- seq(a, b)
        S <- .clone_cp_stats(M[, car, drop = FALSE], lo, hi, cand) * w[car]
        cls$pos[k] <- cand[which.max(colSums(S))]
        cls$spread[k] <- diff(range(cand[apply(S, 1, which.max)]))
    }
    cls
}

# re-call each subcluster on the shared segments: most common direction, then most common state with that direction
.clone_resegment <- function(P, chr, pos) {
    starts <- sort(unique(c(tapply(seq_along(chr), chr, min), pos)))
    seg <- findInterval(seq_along(chr), starts)
    Q <- P
    for (s in colnames(P)) {
        m <- tapply(P[, s], seg, function(x) { d <- .clone_mode(sign(x - 3)); .clone_mode(x[sign(x - 3) == d]) })
        Q[, s] <- m[as.character(seg)]
    }
    list(P = Q, seg = seg)
}

.clone_refine_group <- function(P, M, chr, tol = 30) {
    bk <- .clone_breaks(.clone_dir(P), chr)      # event boundaries: changes of direction
    if (is.null(bk) || !nrow(bk)) {
        return(list(P = P, seg = findInterval(seq_along(chr), tapply(seq_along(chr), chr, min)), breaks = NULL))
    }
    bk <- .clone_cluster_breaks(bk, tol)
    br <- .clone_place_breaks(bk, M, chr)
    rs <- .clone_resegment(P, chr, br$pos)
    list(P = rs$P, seg = rs$seg, breaks = br)
}


## ---- 2. candidate clones ----

.clone_candidates <- function(P2, seg, R, cells_by_sub, min_cells = 5, min_genes = 30, t_min = 5, d_frac = 0.5,
                              min_amp = 0.4, q_clonal = 0.85, max_start = 30) {
    subs <- colnames(P2); ns <- lengths(cells_by_sub)[subs]
    nseg <- table(seg); fs <- as.integer(names(nseg)[nseg >= min_genes]); wg <- as.vector(nseg[as.character(fs)])
    first <- match(fs, seg)
    D <- .clone_dir(P2[first, , drop = FALSE])
    allc <- unlist(cells_by_sub[subs])
    SM <- rowsum(R[seg %in% fs, allc, drop = FALSE], seg[seg %in% fs]) / wg         # segments x cells
    colnames(SM) <- allc
    big <- subs[ns >= min_cells]
    wfrac <- function(L) as.vector(L %*% ns[colnames(L)]) / sum(ns[colnames(L)])
    # baseline per cell: segments neutral in >= q_clonal of the cells
    neut <- wfrac(D[, big, drop = FALSE] == 0) >= q_clonal
    if (any(neut)) SM <- sweep(SM, 2, colSums(wg[neut] * SM[neut, , drop = FALSE]) / sum(wg[neut]))
    # signal amplitude per subcluster (purity of spots, weaker signal) on segments with one shared call
    clon <- wfrac(D[, big, drop = FALSE] == 1) >= q_clonal | wfrac(D[, big, drop = FALSE] == -1) >= q_clonal
    subm <- vapply(subs, function(s) rowMeans(SM[, as.character(cells_by_sub[[s]]), drop = FALSE]), numeric(length(fs)))
    if (!is.matrix(subm)) subm <- matrix(subm, nrow = length(fs), dimnames = list(NULL, subs))
    amp <- setNames(rep(1, length(subs)), subs)
    if (sum(clon) >= 2) { m <- rowMeans(subm[clon, big, drop = FALSE])
        amp <- apply(subm[clon, , drop = FALSE], 2, function(y) sum(wg[clon] * y * m) / sum(wg[clon] * m^2)) }
    low <- subs[amp < min_amp]
    big <- setdiff(big, low)
    if (length(big) < 2) return(list(clone = setNames(rep("A", length(subs)), subs), low = low))
    SM <- sweep(SM, 2, pmax(amp, 0.2)[rep(subs, lengths(cells_by_sub[subs]))], "/")
    subm2 <- vapply(big, function(s) rowMeans(SM[, as.character(cells_by_sub[[s]]), drop = FALSE]), numeric(length(fs)))
    if (!is.matrix(subm2)) subm2 <- matrix(subm2, nrow = length(fs), dimnames = list(NULL, big))
    S1 <- P2[first, big, drop = FALSE]; sc1 <- S1 == 2 | S1 == 4
    step <- if (any(sc1)) .clone_wmedian(abs(subm2[sc1]), outer(wg, ns[big])[sc1]) else 0.2
    d_min <- d_frac * step

    sig <- apply(D[, big, drop = FALSE], 2, paste, collapse = "")
    grp <- setNames(match(sig, unique(sig)), big)
    if (max(grp) > max_start) {
        grp <- setNames(stats::cutree(stats::hclust(stats::dist(t(subm2 * sqrt(wg))), "ward.D2"), k = max_start), big)
    }
    cells_of <- function(ss) as.character(unlist(cells_by_sub[ss]))
    maj <- function(ss) apply(D[, ss, drop = FALSE], 1, .clone_mode)
    distinct <- function(a, b) {
        sa <- names(grp)[grp == a]; sb <- names(grp)[grp == b]
        da <- maj(sa); db <- maj(sb)
        cand <- which(da != db); if (!length(cand)) return(FALSE)
        ca <- cells_of(sa); cb <- cells_of(sb)
        for (i in cand) { x <- SM[i, ca]; y <- SM[i, cb]
            se <- sqrt(stats::var(x) / length(x) + stats::var(y) / length(y)); if (!is.finite(se) || se == 0) se <- 1e-3
            if (abs(mean(x) - mean(y)) >= d_min && abs(mean(x) - mean(y)) / se >= t_min) return(TRUE) }
        FALSE
    }
    ham <- function(a, b) sum(wg * (maj(names(grp)[grp == a]) != maj(names(grp)[grp == b])))
    score <- function(a, b) if (distinct(a, b)) Inf else ham(a, b)
    ids <- sort(unique(grp)); PS <- matrix(Inf, max(ids), max(ids))
    if (length(ids) > 1) { pr <- utils::combn(ids, 2); for (k in seq_len(ncol(pr))) PS[pr[1, k], pr[2, k]] <- score(pr[1, k], pr[2, k]) }
    repeat {                                              # merge the most similar non-distinct pair
        ids <- sort(unique(grp)); if (length(ids) < 2) break
        sp <- PS[ids, ids, drop = FALSE]; if (!any(is.finite(sp))) break
        w <- which(sp == min(sp), arr.ind = TRUE)[1, ]; a <- ids[w[1]]; b <- ids[w[2]]
        grp[grp == b] <- a; PS[b, ] <- Inf; PS[, b] <- Inf
        for (o in setdiff(sort(unique(grp)), a)) PS[min(a, o), max(a, o)] <- score(min(a, o), max(a, o))
    }
    # small and low-signal subclusters join the nearest clone
    cm <- vapply(sort(unique(grp)), function(g) rowMeans(SM[, cells_of(names(grp)[grp == g]), drop = FALSE]), numeric(length(fs)))
    if (!is.matrix(cm)) cm <- matrix(cm, nrow = length(fs))
    gid <- sort(unique(grp))
    for (s in setdiff(subs, big)) { m <- rowMeans(SM[, cells_of(s), drop = FALSE])
        grp[s] <- gid[which.min(colSums(wg * (cm - m)^2))] }
    list(clone = setNames(as.character(grp[subs]), subs), low = low)
}


## ---- 3. parsimony re-call and merging ----

.clone_recall <- function(cons, seg, chr, Mc, w, lambda = 10, n_iter = 3) {
    first <- !duplicated(seg); sg <- seg[first]; ch <- chr[first]
    len <- as.vector(table(seg)[as.character(sg)])
    Y <- rowsum(Mc, seg)[as.character(sg), , drop = FALSE] / len
    S0 <- cons[first, , drop = FALSE]
    neut <- apply(S0, 1, function(x) all(x == 3))
    if (any(neut)) Y <- sweep(Y, 2, colSums(len[neut] * Y[neut, , drop = FALSE]) / sum(len[neut]))
    sig2 <- apply(Mc, 2, function(y) stats::median(tapply(y, seg, stats::var), na.rm = TRUE))
    big <- len >= 30
    cl_ev <- big & apply(S0, 1, function(x) all(x == x[1]) && x[1] != 3)
    amp <- setNames(rep(1, ncol(Y)), colnames(Y))
    if (sum(cl_ev) >= 2) { m <- rowMeans(Y[cl_ev, , drop = FALSE])
        amp <- apply(Y[cl_ev, , drop = FALSE], 2, function(y) max(0.3, sum(len[cl_ev] * y * m) / sum(len[cl_ev] * m^2))) }
    Ya <- sweep(Y, 2, amp, "/")
    mu <- vapply(1:6, function(s) { v <- Ya[big, , drop = FALSE][S0[big, , drop = FALSE] == s]
        if (length(v) >= 2) stats::median(v) else NA_real_ }, numeric(1))
    step <- mean(abs(mu[c(2, 4)]), na.rm = TRUE); if (!is.finite(step)) step <- 0.2
    fb <- c(-2, -1, 0, 1, 2, 4) * step; mu[is.na(mu)] <- fb[is.na(mu)]; mu[3] <- 0
    wk <- w[colnames(Y)]
    b <- rep(0, nrow(Y)); S1 <- S0
    for (it in seq_len(if (ncol(Y) > 1) n_iter else 1)) {
        if (ncol(Y) > 1) b <- as.vector(((Y - sweep(matrix(mu[S1], nrow = nrow(S1)), 2, amp, "*")) %*% wk) / sum(wk))
        Yb <- Y - b; S1 <- S0
        for (c0 in unique(ch)) { j <- which(ch == c0)
            allowed <- lapply(j, function(i) sort(unique(c(3, S0[i, ]))))
            for (k in colnames(S0)) {
                cost <- lapply(seq_along(j), function(t) len[j[t]] * (Yb[j[t], k] - amp[k] * mu[allowed[[t]]])^2 / sig2[k])
                V <- cost[[1]]; bp <- list()
                if (length(j) > 1) for (t in 2:length(j)) {
                    prev <- allowed[[t - 1]]; cur <- allowed[[t]]
                    tr <- outer(seq_along(prev), seq_along(cur), function(a, b) V[a] + lambda * (prev[a] != cur[b]))
                    bp[[t]] <- apply(tr, 2, which.min); V <- apply(tr, 2, min) + cost[[t]]
                }
                path <- integer(length(j)); path[length(j)] <- which.min(V)
                if (length(j) > 1) for (t in length(j):2) path[t - 1] <- bp[[t]][path[t]]
                S1[j, k] <- mapply(function(t, p) allowed[[t]][p], seq_along(j), path)
            }
        }
    }
    list(cons = S1[match(seg, sg), , drop = FALSE], mu = mu, amp = amp)
}

.clone_merge <- function(cl, P2, seg, chr, R, cells_by_sub, lambda = 10, min_event_genes = 100, d_frac = 0.7) {
    clone <- cl$clone; ns <- lengths(cells_by_sub)
    repeat {
        ks <- sort(unique(clone))
        members <- lapply(ks, function(k) names(clone)[clone == k])
        cons <- vapply(members, function(ss) { if (any(!ss %in% cl$low)) ss <- setdiff(ss, cl$low)
            apply(P2[, ss, drop = FALSE], 1, .clone_wmode, w = ns[ss]) }, numeric(nrow(P2)))
        if (!is.matrix(cons)) cons <- matrix(cons, ncol = 1)
        Mc <- vapply(members, function(ss) rowMeans(R[, unlist(cells_by_sub[ss]), drop = FALSE]), numeric(nrow(R)))
        if (!is.matrix(Mc)) Mc <- matrix(Mc, ncol = 1)
        colnames(cons) <- colnames(Mc) <- ks
        w <- setNames(vapply(members, function(ss) sum(ns[ss]), numeric(1)), ks)
        pr <- .clone_recall(cons, seg, chr, Mc, w, lambda = lambda)
        if (length(ks) < 2) break
        D <- .clone_dir(pr$cons); pairs <- utils::combn(ks, 2); step <- mean(abs(pr$mu[c(2, 4)]))
        score <- apply(pairs, 2, function(p) {
            dif <- D[, p[1]] != D[, p[2]]; if (!any(dif)) return(0)
            r <- rle(paste(chr, dif)); e <- cumsum(r$lengths); s <- e - r$lengths + 1
            ev <- which(endsWith(r$values, "TRUE") & r$lengths >= min_event_genes)
            if (!length(ev)) return(0)
            max(vapply(ev, function(i) abs(stats::median(Mc[s[i]:e[i], p[1]] - Mc[s[i]:e[i], p[2]])), numeric(1))) /
                (step * min(pr$amp[p]))
        })
        if (min(score) >= d_frac) break
        p <- pairs[, which.min(score)]
        clone[clone == p[2]] <- p[1]
    }
    sz <- tapply(ns[names(clone)], clone, sum); o <- names(sz)[order(-sz)]
    lab <- setNames(LETTERS[seq_along(o)], o)
    pc <- pr$cons; colnames(pc) <- lab[colnames(pc)]; pc <- pc[, sort(colnames(pc)), drop = FALSE]
    clone <- setNames(lab[clone], names(clone))
    list(clone = clone, pcons = pc, sizes = tapply(ns[names(clone)], clone, sum))
}


## ---- 4. events and the clone tree ----

# events: runs of segments on one chromosome with the same direction and the same carrier clones
.clone_shared_events <- function(cons, seg, chr, min_genes = 30) {
    cl <- colnames(cons); first <- !duplicated(seg); ch <- chr[first]
    len <- as.vector(table(seg)[as.character(seg[first])]); st <- which(first); en <- c(st[-1] - 1L, length(seg))
    rows <- list()
    for (d in c(1, -1)) {
        car <- apply(.clone_dir(cons[first, , drop = FALSE]) == d, 1, function(x) paste(cl[x], collapse = ""))
        r <- rle(paste(ch, car)); e <- cumsum(r$lengths); s <- e - r$lengths + 1
        for (k in seq_along(s)) { if (car[s[k]] == "") next
            rows[[length(rows) + 1]] <- data.frame(chr = ch[s[k]], start = st[s[k]], end = en[e[k]], dir = d,
                                                   genes = sum(len[s[k]:e[k]]), carriers = car[s[k]], stringsAsFactors = FALSE) }
    }
    if (!length(rows)) return(data.frame(chr = character(), start = integer(), end = integer(), dir = numeric(),
                                         genes = integer(), carriers = character(), stringsAsFactors = FALSE))
    ev <- do.call(rbind, rows); ev <- ev[ev$genes >= min_genes, ]
    ev[order(ev$start), ]
}

.clone_conflicts <- function(ev) {
    cs <- unique(ev$carriers); if (length(cs) < 2) return(character(0))
    sets <- strsplit(cs, ""); bad <- character(0)
    for (i in 1:(length(cs) - 1)) for (j in (i + 1):length(cs)) { x <- intersect(sets[[i]], sets[[j]])
        if (length(x) && length(x) < length(sets[[i]]) && length(x) < length(sets[[j]])) bad <- c(bad, cs[i], cs[j]) }
    unique(bad)
}

# tree from nested carrier sets; conflicting carrier sets (fewest genes first) are left off. Branch length = genes
# changed on the branch; leaves get y positions in depth-first order.
.clone_tree_layout <- function(ev, sizes) {
    clones <- sort(names(sizes)); unplaced <- ev[0, ]
    repeat { bad <- .clone_conflicts(ev); if (!length(bad)) break
        g <- tapply(ev$genes[ev$carriers %in% bad], ev$carriers[ev$carriers %in% bad], sum)
        drop <- names(g)[which.min(g)]
        unplaced <- rbind(unplaced, ev[ev$carriers == drop, ]); ev <- ev[ev$carriers != drop, ] }
    sets <- unique(c(paste(clones, collapse = ""), ev$carriers, clones))
    sets <- sets[order(-nchar(sets), sets)]; memb <- strsplit(sets, "")
    parent <- vapply(seq_along(sets), function(i) { if (i == 1) return(NA_character_)
        sup <- which(vapply(seq_along(sets), function(j) j != i && nchar(sets[j]) > nchar(sets[i]) &&
                                all(memb[[i]] %in% memb[[j]]), logical(1)))
        if (!length(sup)) NA_character_ else sets[sup[which.min(nchar(sets[sup]))]] }, character(1))
    lab <- function(e) if (!nrow(e)) "" else paste(sprintf("%s%s:%.0f-%.0fMb", ifelse(e$dir > 0, "+", "-"),
        sub("^chr", "", e$chr), e$start_pos / 1e6, e$end_pos / 1e6), collapse = " ")
    tr <- data.frame(node = sets, parent = parent, stringsAsFactors = FALSE)
    tr$events <- vapply(sets, function(s) lab(ev[ev$carriers == s, ]), character(1))
    tr$len <- vapply(sets, function(s) sum(ev$genes[ev$carriers == s]), numeric(1))
    tr$x <- NA_real_; tr$y <- NA_real_; leaves <- character(0)
    walk <- function(n, x0) { i <- which(tr$node == n); tr$x[i] <<- x0 + tr$len[i]
        ch <- tr$node[which(tr$parent == n)]
        if (!length(ch)) { leaves <<- c(leaves, n); tr$y[i] <<- length(leaves) }
        else { for (c0 in ch[order(-nchar(ch), ch)]) walk(c0, tr$x[i]); tr$y[i] <<- mean(tr$y[match(ch, tr$node)]) } }
    walk(tr$node[1], 0)
    tr$px <- ifelse(is.na(tr$parent), 0, tr$x[match(tr$parent, tr$node)])
    tr$leaf <- !tr$node %in% tr$parent
    tr$size <- ifelse(tr$leaf, as.numeric(sizes[tr$node]), NA_real_)
    tr$truncal <- tr$node == sets[1]
    list(tree = tr, leaves = leaves, unplaced = unplaced)
}


## ---- figure: clone tree + clone copy number ----

.clone_plot <- function(r, chr, HMM_type, title, file) {
    tr <- r$tree; leaves <- tr$node[tr$leaf][order(tr$y[tr$leaf])]
    wrap <- function(s, n = 4) { w <- strsplit(s, " ")[[1]]; if (!length(w)) return("")
        if (length(w) > 12) w <- c(w[1:12], sprintf("(+%d more)", length(w) - 12))
        paste(vapply(split(w, ceiling(seq_along(w) / n)), paste, character(1), collapse = " "), collapse = "\n") }
    tr$lab <- vapply(tr$events, wrap, character(1))
    tr$py <- tr$y[match(tr$parent, tr$node)]
    xmax <- max(tr$x) * 1.3 + 1
    p1 <- ggplot2::ggplot(tr) +
        ggplot2::geom_segment(ggplot2::aes(x = .data$px, xend = .data$x, y = .data$y, yend = .data$y), linewidth = 0.8, colour = "grey30") +
        ggplot2::geom_segment(data = tr[!is.na(tr$parent), ], ggplot2::aes(x = .data$px, xend = .data$px, y = .data$y, yend = .data$py),
                              linewidth = 0.8, colour = "grey30") +
        ggplot2::geom_text(ggplot2::aes(x = .data$px + 0.02 * max(.data$x), y = .data$y, label = .data$lab),
                           hjust = 0, vjust = -0.3, size = 2.4, lineheight = 0.85) +
        ggplot2::geom_point(data = tr[tr$leaf, ], ggplot2::aes(x = .data$x, y = .data$y), size = 2.5) +
        ggplot2::geom_text(data = tr[tr$leaf, ], ggplot2::aes(x = .data$x, y = .data$y, label = sprintf("  clone %s (%d)", .data$node, .data$size)),
                           hjust = 0, size = 3, fontface = "bold") +
        ggplot2::scale_x_continuous(limits = c(0, xmax)) +
        ggplot2::scale_y_reverse(limits = c(length(leaves) + 0.6, 0.4), breaks = NULL) +
        ggplot2::labs(x = "genes changed from normal (trunk: truncal events; branches: derived)", y = NULL) +
        ggplot2::theme_minimal(base_size = 9)
    st <- r$states[, leaves, drop = FALSE]; if (HMM_type == "i3") st <- st + 1
    nc <- 1500; bins <- ceiling(seq_len(nrow(st)) / (nrow(st) / nc))
    mid <- vapply(split(seq_len(nrow(st)), bins), function(i) i[ceiling(length(i) / 2)], 1L)
    df <- data.frame(x = rep(seq_along(mid), length(leaves)), y = rep(seq_along(leaves), each = length(mid)),
                     state = factor(as.vector(st[mid, ]), levels = 1:6))
    cols <- c("1" = "#08306B", "2" = "#6BAED6", "3" = "#F2F2F2", "4" = "#FC9272", "5" = "#DE2D26", "6" = "#67000D")
    labs <- if (HMM_type == "i6") c("0", "1", "2", "3", "4", "≥5") else c("", "loss", "neutral", "gain", "", "")
    cb <- cumsum(table(factor(chr[mid], levels = unique(chr[mid]))))
    cm <- cb - table(factor(chr[mid], levels = unique(chr[mid]))) / 2
    p2 <- ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$state)) + ggplot2::geom_raster() +
        ggplot2::geom_vline(xintercept = cb[-length(cb)] + 0.5, linewidth = 0.2) +
        ggplot2::scale_fill_manual(values = cols, labels = setNames(labs, 1:6), drop = FALSE,
                                   name = if (HMM_type == "i6") "copies" else "") +
        ggplot2::scale_y_reverse(breaks = seq_along(leaves), labels = paste("clone", leaves), expand = c(0, 0)) +
        ggplot2::scale_x_continuous(breaks = as.vector(cm), labels = sub("^chr", "", names(cm)), expand = c(0, 0)) +
        ggplot2::labs(x = NULL, y = NULL, title = "clone copy number") +
        ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(panel.grid = ggplot2::element_blank(),
                                                              legend.position = "bottom")
    grDevices::png(file, width = 15, height = 2 + 0.8 * length(leaves), units = "in", res = 110)
    gridExtra::grid.arrange(p1, p2, ncol = 2, widths = c(1, 1.3),
                            top = sprintf("%s: %d cells in %d clones%s", title, sum(r$sizes), length(r$sizes),
                                          if (nrow(r$unplaced)) "; some events conflict with the tree (see tree file)" else ""))
    grDevices::dev.off()
}
