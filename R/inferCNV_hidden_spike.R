

.build_and_add_hspike <- function(infercnv_obj, sim_method=c('meanvar', 'simple', 'splatter', 'resample'), aggregate_normals=FALSE) {
    
    sim_method = match.arg(sim_method)

    if (sim_method == 'resample') {
        return(.build_and_add_hspike_resample(infercnv_obj, aggregate_normals=aggregate_normals))
    }

    flog.info("Adding h-spike")

    if (has_reference_cells(infercnv_obj)) {
        if (aggregate_normals) {
            ## for experimental use / data exploration
            normal_cells_idx_lists = list();
            normal_cells_idx_lists[[ 'normalsToUse' ]] = unlist(infercnv_obj@reference_grouped_cell_indices)
        } else {
            ## the usual method to use
            normal_cells_idx_lists = infercnv_obj@reference_grouped_cell_indices
        }
    } else {
        ## the reference-less case:
        idx = unlist(infercnv_obj@observation_grouped_cell_indices)
        normal_cells_idx_lists = list()
        normal_cells_idx_lists[[ 'normalsToUse' ]] = idx 
        flog.info("-no normals defined, using all observation cells as proxy") 
        infercnv_obj_tmp = infercnv_obj
        infercnv_obj_tmp@observation_grouped_cell_indices = normal_cells_idx_lists
    }

    params = list()

    ## build a fake genome with fake chromosomes, alternate between 'normal' and 'variable' regions.

    num_cells = 100
    num_genes_per_chr = 400

    num_total_genes = nrow(infercnv_obj@expr.data)

    chr_info <- .get_hspike_chr_info(num_genes_per_chr, num_total_genes)

    gene_order = do.call(rbind, lapply(chr_info, function(x) { data.frame(chr=x$name, start=seq_len(x$ngenes), stop=seq_len(x$ngenes)) }))
    num_genes = nrow(gene_order)
    rownames(gene_order) <- paste0("gene_", seq_len(num_genes))


    genes_means_use_idx = sample(x=seq_len(nrow(infercnv_obj@expr.data)), size=num_genes, replace=TRUE)

    
    ## do for each group of normal cells
    sim.counts.matrix = NULL 
    reference_grouped_cell_indices = list()
    observation_grouped_cell_indices = list()
    cell_counter = 0
    for (normal_type in names(normal_cells_idx_lists)) {

        flog.info(sprintf("-hspike modeling of %s", normal_type))
        
        normal_cells_idx <- normal_cells_idx_lists[[ normal_type ]]
        
        ## sample gene info from the normal data
        normal_cells_expr = infercnv_obj@expr.data[,normal_cells_idx]
        gene_means_orig = rowMeans(normal_cells_expr)
        gene_means = gene_means_orig[genes_means_use_idx]

        #write.table(gene_means, sprintf("gene_means.before.%s",sub(pattern="[^A-Za-z0-9]", replacement="_",  x=normal_type, perl=TRUE)), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE); ## debugging
        
        gene_means[gene_means==0] = 1e-3 # just make small nonzero values

        
        names(gene_means) = rownames(gene_order)        
        
        mean_p0_table <- NULL
        
        ## simulate normals:
        if (sim_method == 'splatter') {
            params = .estimateSingleCellParamsSplatterScrape(counts=infercnv_obj@count.data[,normal_cells_idx])
            params[["nGenes"]] <- num_genes
            params[["nCells"]] <- num_cells
            sim.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(params, use.genes.means=gene_means)
            
            sim_normal_matrix = counts(sim.scExpObj)
        } else if (sim_method == 'simple') {
            ## using simple
            mean_p0_table <- .get_mean_vs_p0_from_matrix(normal_cells_expr)    
            sim_normal_matrix <- .get_simulated_cell_matrix(gene_means,
                                                                       mean_p0_table=mean_p0_table,
                                                                       num_cells=num_cells,
                                                                       common_dispersion=0.1)
        } else if (sim_method == 'meanvar') {
            ## using mean,var trend
            if (exists("infercnv_obj_tmp")) {
                sim_normal_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj_tmp, gene_means, num_cells, include.dropout=TRUE)
            }
            else {
                sim_normal_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj, gene_means, num_cells, include.dropout=TRUE)
            }
        }


        spike_norm_name = sprintf("simnorm_cell_%s", normal_type)
        colnames(sim_normal_matrix) = paste0(spike_norm_name, seq_len(num_cells))
        rownames(sim_normal_matrix) = rownames(gene_order)
        
        ## apply spike-in multiplier vec
        hspike_gene_means = gene_means
        for (info in chr_info) {
            chr_name = info$name
            cnv = info$cnv
            if (cnv != 1) {
                gene_idx = which(gene_order$chr == chr_name)
                hspike_gene_means[gene_idx] =  hspike_gene_means[gene_idx] * cnv
            }
        }
        
        
        if (sim_method == 'splatter') {
            sim_spiked_cnv.scExpObj = .simulateSingleCellCountsMatrixSplatterScrape(params,
                                                                                    use.genes.means=hspike_gene_means)
            sim_spiked_cnv_matrix = counts(sim_spiked_cnv.scExpObj)
        } else if (sim_method == 'simple') {
            ## using simple
            sim_spiked_cnv_matrix <- .get_simulated_cell_matrix(hspike_gene_means,
                                                                           mean_p0_table=mean_p0_table,
                                                                           num_cells=num_cells,
                                                                           common_dispersion=0.1)
        } else if (sim_method == 'meanvar') {
            if (exists("infercnv_obj_tmp")) {
                sim_spiked_cnv_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj_tmp, hspike_gene_means, num_cells, include.dropout=TRUE)
            }
            else {
                sim_spiked_cnv_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj, hspike_gene_means, num_cells, include.dropout=TRUE)
            }
        }

        spike_tumor_name = sprintf("spike_tumor_cell_%s", normal_type)
        
        colnames(sim_spiked_cnv_matrix) = paste0(spike_tumor_name, seq_len(num_cells))
        rownames(sim_spiked_cnv_matrix) = rownames(gene_order)
        
        if (is.null(sim.counts.matrix)) {
            sim.counts.matrix <- cbind(sim_normal_matrix, sim_spiked_cnv_matrix) 
        } else {
            sim.counts.matrix <- cbind(sim.counts.matrix, sim_normal_matrix, sim_spiked_cnv_matrix)
        }
       
        reference_grouped_cell_indices[[ spike_norm_name ]]  = (cell_counter+1):(cell_counter+num_cells)
        cell_counter = cell_counter + num_cells
        
        observation_grouped_cell_indices[[ spike_tumor_name ]] = (cell_counter+1):(cell_counter+num_cells)
        cell_counter = cell_counter + num_cells
    }
    
    .hspike <- new( Class="infercnv",
                   expr.data=sim.counts.matrix,
                   count.data=sim.counts.matrix,
                   gene_order=gene_order,
                   reference_grouped_cell_indices=reference_grouped_cell_indices,
                   observation_grouped_cell_indices=observation_grouped_cell_indices)
    
    validate_infercnv_obj(.hspike)

    .hspike <- normalize_counts_by_seq_depth(.hspike, median(colSums(normal_cells_expr))) # make same target counts/cell as normals

    infercnv_obj@.hspike <- .hspike

    return(infercnv_obj)
}




.get_hspike_chr_info <- function(num_genes_each, num_total) {

    num_remaining = num_total - 10*num_genes_each
    if (num_remaining < num_genes_each) {
        num_remaining = num_genes_each
    }

    ## design for fake chr
    chr_info = list(list(name='chrA',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_0',
                         cnv=0.01,
                         ngenes=num_genes_each),
                    list(name='chr_B',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_0pt5',
                         cnv=0.5,
                         ngenes=num_genes_each),
                    list(name='chr_C',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_1pt5',
                         cnv=1.5,
                         ngenes=num_genes_each),
                    list(name='chr_D',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_2pt0',
                         cnv=2.0,
                         ngenes=num_genes_each),
                    list(name='chr_E',
                         cnv=1,
                         ngenes=num_genes_each),
                    list(name='chr_3pt0',
                         cnv=3,
                         ngenes=num_genes_each),
                    list(name='chr_F',
                         cnv=1,
                         ngenes=num_remaining)
                    )

    return(chr_info)

}


## ---- h-spike from resampled reference cells, at the data's gene density ----
#
# The 'resample' h-spike is built from the real reference cells instead of a fitted mean-variance model:
#  - fake chromosomes have the median number of genes per chromosome of the data (10-400), so smoothing averages over
#    as many genes as it does in the data; each copy-number level is repeated over enough fake chromosomes to give it
#    >= min_genes_per_level genes, and a last neutral chromosome brings the h-spike to the data's number of genes;
#  - each fake chromosome takes a run of consecutive real genes (real local expression structure);
#  - h-spike normal cells are real reference cells; h-spike tumour cells are other real reference cells of the same
#    group whose counts on a fake chromosome of level c are binomially thinned (c < 1) or topped up with the counts
#    of floor(c - 1) other cells plus a thinned (c - 1 - floor(c - 1)) share of one more (c > 1).
# Chromosome names carry the level ('hs<k>_cnv:<c>'), read by .hspike_chr_levels().

.hspike_chr_levels <- function(hspike_obj) {
    chrs <- unique(as.character(hspike_obj@gene_order$chr))
    if (all(grepl("^hs[0-9]+_cnv:", chrs))) {
        return(stats::setNames(as.numeric(sub("^hs[0-9]+_cnv:", "", chrs)), chrs))
    }
    info <- .get_hspike_chr_info(1, 1)
    stats::setNames(vapply(info, function(x) x$cnv, numeric(1)), vapply(info, function(x) x$name, character(1)))
}

.build_and_add_hspike_resample <- function(infercnv_obj, aggregate_normals=FALSE, num_cells=100, min_genes_per_level=400) {

    flog.info("Adding h-spike (resampled reference cells)")

    if (has_reference_cells(infercnv_obj)) {
        normal_cells_idx_lists <- if (aggregate_normals) list(normalsToUse=unlist(infercnv_obj@reference_grouped_cell_indices))
                                  else infercnv_obj@reference_grouped_cell_indices
    } else {
        flog.info("-no normals defined, using all observation cells as proxy")
        normal_cells_idx_lists <- list(normalsToUse=unlist(infercnv_obj@observation_grouped_cell_indices))
    }

    go <- infercnv_obj@gene_order
    real_chr <- as.character(go$chr)
    chr_len <- table(factor(real_chr, levels=unique(real_chr)))
    L <- as.integer(max(10, min(400, round(stats::median(chr_len)))))
    reps <- ceiling(min_genes_per_level / L)
    levels <- c(1, 0.01, 1, 0.5, 1, 1.5, 1, 2, 1, 3)
    num_total <- nrow(infercnv_obj@expr.data)
    lens <- rep(L, reps * length(levels))
    lev <- rep(levels, reps)
    rest <- max(L, num_total - sum(lens))
    lens <- c(lens, rest); lev <- c(lev, 1)
    chr_names <- sprintf("hs%03d_cnv:%g", seq_along(lens), lev)
    flog.info(sprintf("-h-spike: %d fake chromosomes of %d genes (data median %g genes per chromosome), %d per level",
                      length(lens) - 1, L, stats::median(chr_len), reps))

    ## runs of consecutive real genes, one per fake chromosome
    chr_start <- cumsum(c(0, as.vector(chr_len)))[seq_along(chr_len)]
    gene_idx <- unlist(lapply(lens, function(len) {
        ok <- which(chr_len >= len)
        if (length(ok)) {
            k <- if (length(ok) == 1) ok else sample(ok, 1)
            chr_start[k] + sample.int(chr_len[k] - len + 1, 1) - 1 + seq_len(len)
        } else {
            sample.int(num_total, len, replace=TRUE)
        }
    }))
    lev_gene <- rep(lev, lens)
    gene_order <- data.frame(chr=rep(chr_names, lens), start=unlist(lapply(lens, seq_len)), stop=unlist(lapply(lens, seq_len)))
    rownames(gene_order) <- paste0("gene_", seq_len(nrow(gene_order)))

    counts <- infercnv_obj@count.data
    thin <- function(x, p) { x <- round(as.matrix(x)); matrix(stats::rbinom(length(x), x, p), nrow=nrow(x)) }

    sim.counts.matrix <- NULL
    reference_grouped_cell_indices <- list()
    observation_grouped_cell_indices <- list()
    cell_counter <- 0
    for (normal_type in names(normal_cells_idx_lists)) {
        flog.info(sprintf("-hspike modeling of %s", normal_type))
        idx <- normal_cells_idx_lists[[ normal_type ]]
        pick <- function(k) if (length(idx) >= k) sample(idx, k) else sample(idx, k, replace=TRUE)
        two <- pick(2 * num_cells)
        nrm <- two[seq_len(num_cells)]; tum <- two[num_cells + seq_len(num_cells)]
        sim_normal_matrix <- as.matrix(counts[gene_idx, nrm, drop=FALSE])
        sim_tumor_matrix <- as.matrix(counts[gene_idx, tum, drop=FALSE])
        for (c0 in setdiff(unique(lev_gene), 1)) {
            g <- which(lev_gene == c0)
            if (c0 < 1) {
                sim_tumor_matrix[g, ] <- thin(sim_tumor_matrix[g, , drop=FALSE], c0)
            } else {
                extra <- c0 - 1
                for (w in seq_len(floor(extra))) sim_tumor_matrix[g, ] <- sim_tumor_matrix[g, ] + as.matrix(counts[gene_idx[g], pick(num_cells), drop=FALSE])
                fr <- extra - floor(extra)
                if (fr > 0) sim_tumor_matrix[g, ] <- sim_tumor_matrix[g, ] + thin(counts[gene_idx[g], pick(num_cells), drop=FALSE], fr)
            }
        }
        spike_norm_name <- sprintf("simnorm_cell_%s", normal_type)
        spike_tumor_name <- sprintf("spike_tumor_cell_%s", normal_type)
        colnames(sim_normal_matrix) <- paste0(spike_norm_name, seq_len(num_cells))
        colnames(sim_tumor_matrix) <- paste0(spike_tumor_name, seq_len(num_cells))
        rownames(sim_normal_matrix) <- rownames(sim_tumor_matrix) <- rownames(gene_order)
        sim.counts.matrix <- cbind(sim.counts.matrix, sim_normal_matrix, sim_tumor_matrix)
        reference_grouped_cell_indices[[ spike_norm_name ]] <- (cell_counter + 1):(cell_counter + num_cells)
        cell_counter <- cell_counter + num_cells
        observation_grouped_cell_indices[[ spike_tumor_name ]] <- (cell_counter + 1):(cell_counter + num_cells)
        cell_counter <- cell_counter + num_cells
    }

    .hspike <- new(Class="infercnv",
                   expr.data=sim.counts.matrix,
                   count.data=sim.counts.matrix,
                   gene_order=gene_order,
                   reference_grouped_cell_indices=reference_grouped_cell_indices,
                   observation_grouped_cell_indices=observation_grouped_cell_indices)
    validate_infercnv_obj(.hspike)

    normal_cells_expr <- infercnv_obj@expr.data[, unlist(normal_cells_idx_lists), drop=FALSE]
    .hspike <- normalize_counts_by_seq_depth(.hspike, median(colSums(normal_cells_expr)))
    infercnv_obj@.hspike <- .hspike
    return(infercnv_obj)
}


##
.sim_foreground <- function(infercnv_obj, sim_method) {

    flog.info("## simulating foreground")

    expr.matrix <- infercnv_obj@expr.data

    samples <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)


    normal_data <- expr.matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)]
    target_total_counts <- median(colSums(normal_data))

    params <- NULL
    if (sim_method == 'splatter') {
        params <- .estimateSingleCellParamsSplatterScrape(infercnv_obj@count.data[, unlist(infercnv_obj@reference_grouped_cell_indices)])
    }


    mean_p0_table <- NULL
    if (sim_method == 'simple') {
        mean_p0_table <- .get_mean_vs_p0_from_matrix(normal_data)
    }
    
    for (sample_name in names(samples)) {

        cell_idx = samples[[ sample_name ]]

        sample_expr_data = expr.matrix[, cell_idx]
        gene_means = rowMeans(sample_expr_data)
        gene_means[gene_means==0] <- 1e-3 #avoid zeros, breaks splatter sim

        num_cells = ncol(sample_expr_data)

        ## sim the tumor matrix
        if (sim_method == 'simple') {
            
            sim_matrix <- .get_simulated_cell_matrix(gene_means,
                                                                mean_p0_table=mean_p0_table,
                                                                num_cells=num_cells,
                                                                common_dispersion=0.1)
        } else if (sim_method == 'splatter') {
            params[['nCells']] <- num_cells
            sim_matrix <- .simulateSingleCellCountsMatrixSplatterScrape(params, gene_means)
            sim_matrix <- counts(sim_matrix)
        } else if (sim_method == 'meanvar') {
            sim_matrix <- .get_simulated_cell_matrix_using_meanvar_trend(infercnv_obj,
                                                                         gene_means,
                                                                         num_cells,
                                                                         include.dropout=TRUE)
        } else {
            stop(sprintf("not recognizing --sim_method: %s", args$sim_method))
        }

        expr.matrix[, cell_idx] <- sim_matrix

    }

    infercnv_obj@expr.data <- expr.matrix

    infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj, target_total_counts)

    return(infercnv_obj)
}

