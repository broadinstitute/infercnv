
#' @title infercnv_auto_reference
#'
#' @description Finds normal (copy-number neutral) cells to use as the reference when none are known, and classifies
#' the remaining cells as malignant or normal. Two passes of inferCNV up to the subclustering step (no HMM):
#'
#' 1. Provisional reference. All cells are subclustered together (Leiden; no reference). The provisional reference
#'    is taken from (reference_method = "auto": the first that applies):
#'    (a) normal_groups, annotation groups known to be normal;
#'    (b) "markers": subclusters expressing markers of non-malignant lineages (immune, endothelial, fibroblast;
#'        normal_markers), as SCEVAN does to find confident normal cells (De Falco et al. 2023, Nat Commun 14:1074):
#'        per panel, the subclusters are split in two by mean log panel expression (k-means), and the high group is
#'        marker-positive when it is at least marker_fold higher;
#'    (c) "self_centered": the expression cluster with the flattest self-centred profile, the principle of CopyKAT's
#'        baseline search (Gao et al. 2021, Nat Biotechnol 39:599). Cells are first clustered by expression (standard
#'        Seurat workflow: log normalisation, flat_n_features variable genes, PCA, shared nearest neighbours, Louvain at
#'        flat_resolution), so that clusters correspond to cell types. Each cluster's mean log2(CPM + 1) profile is
#'        centred on its own genome-wide mean (not on the mean of all cells, which a tumour majority would dominate),
#'        smoothed along each chromosome (running mean of window_length genes) and scored by its mean square. The
#'        flattest cluster is the provisional reference; clusters scoring within flat_tolerance times its score are
#'        added (flat_tolerance = 1: the flattest only; larger values admit more clusters for exploring which cell
#'        types look normal). All scores are written to auto_reference.flatness.tsv. The score is not CNV-specific
#'        (each cell type's expression pattern along the genome contributes), so it is reliable for tumours with many
#'        large CNVs; in tumours with few CNVs (e.g. some gliomas) a tumour cluster can be the flattest. Check the
#'        table.
#' 2. inferCNV is rerun with the provisional reference. Each subcluster is scored as in Tirosh et al. 2016 (Science
#'    352:189) and Puram et al. 2017 (Cell 171:1611): CNV score (mean squared deviation of its mean profile) and CNV
#'    correlation (Pearson correlation of its mean profile with the mean profile of the top-scoring subclusters,
#'    together at least 5% of the cells). Thresholds come from random groups of reference cells of the same size
#'    (score_quantile of their scores and correlations). Subclusters above both thresholds are malignant, below both
#'    normal, otherwise unresolved. Classification is per subcluster, as single cells are too noisy (see
#'    HMM_min_subcluster_size in run()). If most cells called malignant are marker-positive, the reference is probably
#'    inverted and a warning is given.
#'
#' Marker panels assume the tumour itself is not of an immune, endothelial or fibroblast lineage (e.g. not lymphoma or
#' sarcoma); give normal_groups, normal_markers, or normal_markers = NULL otherwise.
#'
#' The returned object has the provisional reference cells and the subclusters classified as normal as reference
#' groups (named "normal:<original group>"), and the other cells as observations under their original groups.
#'
#' @param infercnv_obj infercnv object (references, if any, are ignored and treated as observations)
#' @param out_dir directory for the two passes (subdirectories auto_reference_pass1 and _pass2)
#' @param cutoff expression cutoff passed to run() (default "auto"; see run())
#' @param normal_groups annotation groups known to be normal (optional; used as the provisional reference)
#' @param normal_markers named list of marker gene panels of non-malignant lineages (default: immune, endothelial,
#'                       fibroblast), or NULL to skip the marker step
#' @param marker_fold minimum fold difference of mean panel expression between marker-positive and other subclusters
#'                    (default: 4)
#' @param score_quantile quantile of the reference-group null used as the threshold (default: 0.99)
#' @param min_ref_cells minimum number of provisional reference cells (default: 20)
#' @param reference_method how the provisional reference is found when normal_groups is not given: "auto" (markers
#'                         if enough marker-positive cells, otherwise self_centered), "markers" or "self_centered"
#' @param flat_tolerance self_centered: subclusters whose flatness score is at most flat_tolerance times that of the
#'                       flattest subcluster join the provisional reference (default: 1, the flattest only)
#' @param min_flat_cells self_centered: minimum cells for a cluster to be a candidate (default: 10)
#' @param flat_resolution self_centered: clustering resolution of the expression clusters (default: 1)
#' @param flat_n_features self_centered: number of variable genes for the expression clustering (default: 2000)
#' @param ... further arguments to run() for both passes (e.g. num_threads, k_nn, leiden_resolution)
#'
#' @return list(infercnv_obj = object with the new reference groups, cells = data.frame per cell: cell, group,
#'         subcluster, class (reference / normal / malignant / unresolved), cnv_score, cnv_cor, marker_positive,
#'         pass1_subcluster, expr_cluster), reference_method ("normal_groups", "markers" or "self_centered"), flatness =
#'         data.frame per expression cluster: cluster, cells, flat_score, ratio_to_flattest, marker_positive, selected)
#'
#' @export

infercnv_auto_reference <- function(infercnv_obj, out_dir, cutoff = "auto", normal_groups = NULL,
                                    normal_markers = .default_normal_markers, marker_fold = 4, score_quantile = 0.99,
                                    min_ref_cells = 20, reference_method = c("auto", "markers", "self_centered"),
                                    flat_tolerance = 1, min_flat_cells = 10, flat_resolution = 1, flat_n_features = 2000,
                                    ...) {
    reference_method <- match.arg(reference_method)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    groups <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    cells <- colnames(infercnv_obj@expr.data)
    ann <- data.frame(group = rep(names(groups), lengths(groups)), row.names = cells[unlist(groups)],
                      stringsAsFactors = FALSE)
    counts <- infercnv_obj@count.data
    go <- infercnv_obj@gene_order
    make_obj <- function(ref_cells) {
        a <- ann
        if (length(ref_cells)) a[ref_cells, "group"] <- paste0("normal:", a[ref_cells, "group"])
        refs <- unique(a[ref_cells, "group"])
        CreateInfercnvObject(raw_counts_matrix = counts, annotations_file = a, gene_order_file = go,
                             ref_group_names = if (length(refs)) refs else NULL,
                             chr_exclude = character(0), min_max_counts_per_cell = c(1, Inf))
    }
    upto <- function(obj, sub_dir, by_groups) {
        infercnv::run(obj, cutoff = cutoff, out_dir = file.path(out_dir, sub_dir), cluster_by_groups = by_groups,
                      analysis_mode = "subclusters", HMM = FALSE, denoise = FALSE, no_plot = TRUE, no_prelim_plot = TRUE,
                      save_rds = FALSE, write_expr_matrix = FALSE, up_to_step = 17, ...)
    }
    score <- function(P) colMeans((P - 1)^2)

    ## pass 1: no reference; subclusters of all cells
    flog.info("auto reference, pass 1: no reference, all cells subclustered together")
    o1 <- upto(make_obj(character(0)), "auto_reference_pass1", FALSE)
    sub1 <- unlist(o1@tumor_subclusters$subclusters, recursive = FALSE)
    cn1 <- colnames(o1@expr.data)
    mpos <- .marker_positive_subclusters(counts, sub1, cn1, normal_markers, marker_fold)
    win <- if (!is.null(list(...)$window_length)) list(...)$window_length else 101
    ecl <- .expression_clusters(counts[, cn1, drop = FALSE], flat_resolution, flat_n_features)
    esub <- split(seq_along(cn1), ecl)
    fs <- .self_centered_flatness(counts, o1@gene_order, esub, cn1, win)
    ok <- lengths(esub) >= min_flat_cells
    if (!any(ok)) stop(sprintf("auto reference: no expression cluster of at least %d cells", min_flat_cells))
    fmin <- min(fs[ok])
    flat_sel <- ok & fs <= flat_tolerance * fmin
    flatness <- data.frame(cluster = names(esub), cells = lengths(esub), flat_score = fs, ratio_to_flattest = fs / fmin,
                           marker_positive = .marker_positive_subclusters(counts, esub, cn1, normal_markers, marker_fold),
                           selected = flat_sel, stringsAsFactors = FALSE)
    flatness <- flatness[order(flatness$flat_score), ]
    utils::write.table(flatness, file.path(out_dir, "auto_reference.flatness.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    use_markers <- reference_method == "markers" || (reference_method == "auto" && sum(lengths(sub1)[mpos]) >= min_ref_cells)
    if (!is.null(normal_groups)) {
        method <- "normal_groups"
        ref_cells <- rownames(ann)[ann$group %in% normal_groups]
        flog.info(sprintf("auto reference, pass 1: %d cells of the normal groups %s", length(ref_cells),
                          paste(normal_groups, collapse = ", ")))
    } else if (use_markers) {
        method <- "markers"
        ref_cells <- cn1[unlist(sub1[mpos])]
        flog.info(sprintf("auto reference, pass 1: %d of %d subclusters (%d cells) express non-malignant lineage markers",
                          sum(mpos), length(sub1), length(ref_cells)))
    } else {
        method <- "self_centered"
        ref_cells <- cn1[unlist(esub[flat_sel])]
        flog.info(sprintf(paste0("auto reference, pass 1: %d of %d expression clusters (%d cells) with the flattest self-centred ",
                                 "profiles (flat_tolerance %.2f; next flattest at %.2f x). Not CNV-specific: check ",
                                 "auto_reference.flatness.tsv, or give normal_groups."),
                          sum(flat_sel), length(esub), length(ref_cells), flat_tolerance,
                          if (sum(ok & !flat_sel)) min(fs[ok & !flat_sel]) / fmin else NA))
    }
    if (length(ref_cells) < min_ref_cells) {
        stop(sprintf("auto reference: only %d provisional reference cells; provide reference cells", length(ref_cells)))
    }

    ## pass 2: provisional reference; classify subclusters against random reference groups
    flog.info("auto reference, pass 2: provisional reference")
    o2 <- upto(make_obj(ref_cells), "auto_reference_pass2", TRUE)
    X <- o2@expr.data
    rg <- o2@reference_grouped_cell_indices
    ref_idx <- unlist(rg)
    sub2 <- unlist(o2@tumor_subclusters$subclusters, recursive = FALSE)
    obs2 <- sub2[!vapply(sub2, function(ix) all(ix %in% ref_idx), logical(1))]
    P2 <- vapply(obs2, function(ix) rowMeans(X[, ix, drop = FALSE]), numeric(nrow(X)))
    s2 <- score(P2)
    n2 <- lengths(obs2)
    top <- order(s2, decreasing = TRUE)
    top <- top[seq_len(which(cumsum(n2[top]) >= 0.05 * sum(n2))[1])]
    tprof <- rowMeans(P2[, top, drop = FALSE])
    r2 <- as.numeric(stats::cor(P2, tprof))
    if (exists(".Random.seed", envir = globalenv())) {
        old_seed <- get(".Random.seed", envir = globalenv())
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
    }
    set.seed(1)
    null_at <- function(n) {
        n <- min(n, max(lengths(rg)) %/% 2)
        ok <- which(lengths(rg) >= 2 * n)
        v <- t(vapply(seq_len(100), function(i) {
            g <- ok[sample.int(length(ok), 1, prob = lengths(rg)[ok])]
            p <- rowMeans(X[, rg[[g]][sample.int(length(rg[[g]]), n)], drop = FALSE])
            c(mean((p - 1)^2), stats::cor(p, tprof))
        }, numeric(2)))
        apply(v, 2, stats::quantile, probs = score_quantile, na.rm = TRUE)
    }
    nulls <- vapply(n2, null_at, numeric(2))
    cls <- ifelse(s2 > nulls[1, ] & r2 > nulls[2, ], "malignant",
                  ifelse(s2 <= nulls[1, ] & r2 <= nulls[2, ], "normal", "unresolved"))
    mpos2 <- .marker_positive_subclusters(counts, obs2, colnames(X), normal_markers, marker_fold)
    if (method != "markers" && any(cls == "malignant") && sum(n2[cls == "malignant" & mpos2]) > 0.5 * sum(n2[cls == "malignant"])) {
        flog.warn(paste0("auto reference: most cells called malignant express non-malignant lineage markers; the reference ",
                         "is probably tumour (inverted profiles). Give normal_groups or ref_group_names."))
    }
    flog.info(sprintf("auto reference, pass 2: observation subclusters: %s",
                      paste(sprintf("%s %d (%d cells)", names(table(cls)), as.integer(table(cls)),
                                    as.integer(tapply(n2, cls, sum))), collapse = ", ")))

    cn <- colnames(X)
    res <- data.frame(cell = cn, group = ann[cn, "group"], subcluster = NA_character_, class = NA_character_,
                      cnv_score = NA_real_, cnv_cor = NA_real_, marker_positive = NA, pass1_subcluster = NA_character_,
                      expr_cluster = NA_character_, stringsAsFactors = FALSE)
    rownames(res) <- cn
    res[cn1[unlist(sub1)], "pass1_subcluster"] <- rep(names(sub1), lengths(sub1))
    res[cn1, "expr_cluster"] <- as.character(ecl)
    res[cn[ref_idx], "class"] <- "reference"
    for (j in seq_along(obs2)) {
        cc <- cn[obs2[[j]]]
        res[cc, "subcluster"] <- names(obs2)[j]; res[cc, "class"] <- cls[j]
        res[cc, "cnv_score"] <- s2[j]; res[cc, "cnv_cor"] <- r2[j]; res[cc, "marker_positive"] <- mpos2[j]
    }
    utils::write.table(res, file.path(out_dir, "auto_reference.cells.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    new_ref <- res$cell[res$class %in% c("reference", "normal")]
    list(infercnv_obj = make_obj(new_ref), cells = res, reference_method = method, flatness = flatness)
}


# Expression clusters (standard Seurat workflow; Louvain on the shared-nearest-neighbour graph): cell types, not
# copy-number groups, so that a cluster's self-centred profile is not a mixture of cell types.
.expression_clusters <- function(counts, resolution, n_features) {
    so <- Seurat::CreateSeuratObject(counts, min.cells = 0, min.features = 0)
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so <- Seurat::FindVariableFeatures(so, nfeatures = n_features, verbose = FALSE)
    so <- Seurat::ScaleData(so, verbose = FALSE)
    npcs <- min(30, ncol(so) - 1)
    so <- Seurat::RunPCA(so, npcs = npcs, verbose = FALSE)
    so <- Seurat::FindNeighbors(so, dims = seq_len(npcs), verbose = FALSE)
    so <- Seurat::FindClusters(so, resolution = resolution, verbose = FALSE)
    cl <- as.character(Seurat::Idents(so))
    stats::setNames(paste0("expr_", cl), colnames(so))[colnames(counts)]
}


# Self-centred flatness per subcluster (CopyKAT's baseline principle, no shared baseline): mean log2(CPM + 1) over the
# subcluster's cells (genes of the processed object, in genome order), minus its own genome-wide mean, smoothed by a
# running mean of `window` genes within each chromosome (shorter at chromosome ends); score = mean square.
.self_centered_flatness <- function(counts, gene_order, subs, cell_names, window) {
    g <- intersect(rownames(gene_order), rownames(counts))
    C <- counts[g, cell_names, drop = FALSE]
    N <- methods::as(C %*% Matrix::Diagonal(x = 1e6 / pmax(Matrix::colSums(C), 1)), "CsparseMatrix")
    N@x <- log2(1 + N@x)
    chr <- as.character(gene_order[g, "chr"])
    h <- window %/% 2
    run_mean <- function(x) {
        n <- length(x); cs <- c(0, cumsum(x)); lo <- pmax(1, seq_len(n) - h); hi <- pmin(n, seq_len(n) + h)
        (cs[hi + 1] - cs[lo]) / (hi - lo + 1)
    }
    vapply(subs, function(ix) {
        p <- Matrix::rowMeans(N[, ix, drop = FALSE]); p <- p - mean(p)
        sm <- unlist(lapply(split(p, factor(chr, levels = unique(chr))), run_mean))
        mean(sm^2)
    }, numeric(1))
}


.default_normal_markers <- list(
    immune = c("PTPRC", "CD3E", "CD3D", "CD2", "CD79A", "MS4A1", "LYZ", "CD68", "CD14", "CD163", "C1QA", "C1QB",
               "TYROBP", "FCER1G", "NKG7"),
    endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "KDR", "ESAM"),
    fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "FBLN1"))

# marker-positive subclusters: per panel, mean over a subcluster's cells of log1p(panel counts per 10k); subclusters
# split in two by k-means on this score; the high group is positive if its centre is at least marker_fold higher
.marker_positive_subclusters <- function(counts, subs, cell_names, markers, marker_fold) {
    pos <- rep(FALSE, length(subs))
    if (is.null(markers) || length(subs) < 2) return(pos)
    lib <- Matrix::colSums(counts)
    for (panel in names(markers)) {
        g <- intersect(markers[[panel]], rownames(counts))
        if (length(g) < 2) next
        ps <- log1p(Matrix::colSums(counts[g, , drop = FALSE]) / pmax(lib, 1) * 1e4)
        sc <- vapply(subs, function(ix) mean(ps[cell_names[ix]]), numeric(1))
        if (length(unique(sc)) < 2) next
        km <- stats::kmeans(sc, centers = range(sc))
        hi <- which.max(km$centers)
        if (expm1(max(km$centers)) >= marker_fold * expm1(min(km$centers)) + 1e-9 && max(km$centers) - min(km$centers) >= log(marker_fold)) {
            pos <- pos | km$cluster == hi
        }
    }
    pos
}
