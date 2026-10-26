
#' @title get_spike_dists
#'
#' @description determines the N(mean,sd) parameters for each of the CNV states based on
#' the in silico spike in data (hspike).
#'
#' @param hspike_obj hidden spike object
#'
#' @return cnv_mean_sd list
#'
#' @keywords internal
#' @noRd
#'

get_spike_dists <- function(hspike_obj) {
    
    if (is.null(hspike_obj)) {
        flog.error("get_spike_dists(hspike_obj): Error, hspike obj is null")
        stop("Error")
    }
    if (! is.null(hspike_obj@.hspike)) {
        flog.error("get_spike_dists() should have an hspike obj as param, and this doesn't look like an hspike")
    }
    
    gene_expr_by_cnv = .get_gene_expr_by_cnv(hspike_obj)
        
    cnv_mean_sd = .get_gene_expr_mean_sd_by_cnv(gene_expr_by_cnv)

    return(cnv_mean_sd)
    
}


#' @title .get_gene_expr_by_cnv
#'
#' @description  builds a list containing all intensities corresponding to each of the
#' spiked-in cnv levels.
#'
#' @param hspike_obj
#'
#' @return gene_expr_by_cnv  list, keyed on cnv level, value as vector of residual expr intensities
#'
#' @noRd

.get_gene_expr_by_cnv <- function(hspike_obj) {
    
    chr_levels = .hspike_chr_levels(hspike_obj)

    spike_cell_idx = unlist(hspike_obj@observation_grouped_cell_indices)

    spike.expr.data = hspike_obj@expr.data[, spike_cell_idx]

    gene_expr_by_cnv = list()
    for (chr_name in names(chr_levels)) {
        cnv = sprintf("cnv:%g", chr_levels[[chr_name]])
        chr_gene_idx = which(hspike_obj@gene_order$chr == chr_name)
        gene.expr = c(spike.expr.data[chr_gene_idx, ])
        
        if (cnv %in% names(gene_expr_by_cnv)) {
            gene_expr_by_cnv[[cnv]] = c(gene_expr_by_cnv[[cnv]], gene.expr)
        }  else {
            gene_expr_by_cnv[[cnv]] = gene.expr
        }
    }
    
    return(gene_expr_by_cnv)
}




#' @title .get_gene_expr_mean_sd_by_cnv
#'
#' @description extracts the N(mean,sd) params for each of the cnv levels based
#' on the list of residual expr intensities for each cnv level
#'
#' @param gene_expr_by_cnv data table
#'
#' @return cnv_mean_sd list
#'
#' @noRd

.get_gene_expr_mean_sd_by_cnv <- function(gene_expr_by_cnv) {

    cnv_mean_sd = list()
    
    for (cnv_level in names(gene_expr_by_cnv) ) {
        gene_expr = gene_expr_by_cnv[[ cnv_level ]]

        gene_expr_mean = mean(gene_expr)
        gene_expr_sd = sd(gene_expr)

        cnv_mean_sd[[ cnv_level ]] = list(mean=gene_expr_mean, sd=gene_expr_sd)
    }
    
    return(cnv_mean_sd)
    
}


#' @title .plot_gene_expr_by_cnv
#'
#' @description generates a plot showing the residual expresion intensity distribution along with the
#' reference theoretical densities for each of the corresponding parameterized normal distributions.
#'
#' @param gene_expr_by_cnv list
#'
#' @param cnv_mean_sd list
#'
#' @return ggplot2_plot
#'
#' @noRd

.plot_gene_expr_by_cnv <- function(gene_expr_by_cnv, cnv_mean_sd) {
    
    df  = do.call(rbind, lapply(names(gene_expr_by_cnv), function(x) { data.frame(cnv=x, expr=gene_expr_by_cnv[[x]]) }))

    p = df %>% ggplot(aes(expr,  fill='cnv', colour='cnv'))  +  geom_density(alpha=0.1)

    p = p +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd)) +
        stat_function(fun=dnorm, color='black', args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd)) 

    return(p)
}


#' @title get_hspike_cnv_mean_sd_trend_by_num_cells_fit
#'
#' @description determine the number of cells - to - variance fit for each of the cnv levels.
#'
#' Different numbers of cells are randomly selected from the distribution of residual intensitites at each
#' corresponding CNV level, the variance is computed, and a linear model is then fit.
#'
#' Note, this is similar to what is done in HoneyBadger, but has many differences to how they're doing it there,
#' which appears to involve cnv block length rather than cell number.  Here, block size is not relevant, but rather
#' the number of cells in a pre-defined tumor subcluster.  Also, values are extracted from our in silico spike-in.
#'
#' @param hspike_obj hidden spike object
#'
#' @param plot (boolean flag, default FALSE)
#'
#' @return cnv_level_to_mean_sd_fit list
#'
#' @keywords internal
#' @noRd
#'

get_hspike_cnv_mean_sd_trend_by_num_cells_fit <- function(hspike_obj, plot=FALSE) {
    
    if (all(grepl("^hs[0-9]+_cnv:", unique(as.character(hspike_obj@gene_order$chr))))) {
        return(.hspike_whole_cell_sd_fit(hspike_obj))
    }

    gene_expr_by_cnv <- .get_gene_expr_by_cnv(hspike_obj)
    cnv_level_to_mean_sd = list()
    
    for (cnv_level in names(gene_expr_by_cnv) ) {
        expr_vals = gene_expr_by_cnv[[ cnv_level ]]
        nrounds = 100

        sds <- vapply(seq_len(100), function(ncells) {
            vals <- replicate(nrounds, sample(expr_vals, size=ncells, replace=TRUE))
            if ("matrix" %in% is(vals)) {
                means <- Matrix::rowMeans(vals)
            }
            else {
                means <- mean(vals)
            }
            sd(means)
        }, numeric(1))

        cnv_level_to_mean_sd[[ cnv_level ]] <- sds
    }
    
    if (plot) {
        df = do.call(rbind, lapply(names(cnv_level_to_mean_sd), function(cnv_level) {
            sd=cnv_level_to_mean_sd[[ cnv_level ]]
            data.frame(cnv=cnv_level, num_cells=seq_along(sd), sd=sd)
        }))

        p  = df %>% ggplot(aes_string(x=log(df["num_cells"]),y=log(df["sd"]), color=df["cnv"])) + geom_point()
        pdf("hspike_obs_mean_var_trend.pdf")
        flog.info("plotting: hspike_obs_mean_var_trend.pdf")
        plot(p)
        dev.off()
    }

    ## fit linear model
    # cnv_level_to_mean_sd_fit = list()
    # for (cnv_level in names(cnv_level_to_mean_sd) ) {
    #     sd_vals=cnv_level_to_mean_sd[[ cnv_level ]]
    #     num_cells = seq_along(sd_vals)

    #     flog.info(sprintf("fitting num cells vs. variance for cnv level: %s", cnv_level))
    #     fit = lm(log(sd_vals) ~ log(num_cells)) #note, hbadger does something similar, but not for the hmm cnv state levels
    #     cnv_level_to_mean_sd_fit[[ cnv_level ]] = fit
    # }

    tmp_names <- names(cnv_level_to_mean_sd)
    cnv_level_to_mean_sd_fit <- lapply(tmp_names, function(cnv_level) {
        sd_vals=cnv_level_to_mean_sd[[ cnv_level ]]
        num_cells = seq_along(sd_vals)
        fit = lm(log(sd_vals) ~ log(num_cells))
        fit
    })
    names(cnv_level_to_mean_sd_fit) <- tmp_names
    
    return(cnv_level_to_mean_sd_fit)
    
}




# sd of an n-cell mean at one level: lm of log(sd) on log(num_cells) (i.i.d. resampling of h-spike values), or the
# whole-cell model below (sd^2 = a / n + b)
.hspike_sd_at <- function(fit, num_cells) {
    if (inherits(fit, "hspike_sd_model")) {
        return(sqrt(fit$a / num_cells + fit$b))
    }
    exp(predict(fit, newdata=data.frame(num_cells=num_cells))[[1]])
}

# For the 'resample' h-spike: per level, the spread across genes of the mean profile of n h-spike tumour cells drawn
# as whole cells (keeps per-gene structure and the correlation between neighbouring smoothed genes), for n = 1-100,
# fitted as sd^2 = a / n + b: a sampling part that shrinks with more cells and a part that does not. The 100-cell pool's
# own sampling noise (a / pool size) is taken out of b.
.hspike_whole_cell_sd_fit <- function(hspike_obj, n_max=100, nrounds=30) {
    chr_levels <- .hspike_chr_levels(hspike_obj)
    cells <- unlist(hspike_obj@observation_grouped_cell_indices)
    E <- hspike_obj@expr.data[, cells, drop=FALSE]
    lev_gene <- chr_levels[as.character(hspike_obj@gene_order$chr)]
    fits <- lapply(sort(unique(chr_levels)), function(lv) {
        X <- E[which(lev_gene == lv), , drop=FALSE]
        mu <- mean(X)
        ns <- unique(round(exp(seq(0, log(n_max), length.out=25))))
        v <- vapply(ns, function(n) mean(replicate(nrounds, {
            m <- rowMeans(X[, sample(ncol(X), n, replace=TRUE), drop=FALSE]); mean((m - mu)^2) })), numeric(1))
        f <- stats::lm(v ~ I(1 / ns))
        a <- unname(stats::coef(f)[2]); b <- unname(stats::coef(f)[1])
        if (!is.finite(a) || a < 0) a <- 0
        b <- max(0, b - a / ncol(X))
        structure(list(a=a, b=b, n=ns, var=v), class="hspike_sd_model")
    })
    names(fits) <- sprintf("cnv:%g", sort(unique(chr_levels)))
    fits
}


# The HMM works on expression smoothed over window_length genes. When most chromosomes have fewer genes than one
# window (sparse data: few genes pass the expression cutoff), the smoothed profile of a subcluster wanders by more than
# the gap between the HMM states and the HMM calls that noise as CNVs. Stop (or warn) when the median chromosome has
# fewer than window_length genes; warn below 2 x window_length; list chromosomes with fewer than window_length / 2.
.check_hmm_gene_density <- function(infercnv_obj, window_length, action=c("stop", "warn")) {
    action <- match.arg(action)
    chr <- as.character(infercnv_obj@gene_order$chr)
    n <- table(factor(chr, levels=unique(chr)))
    med <- stats::median(n)
    msg <- sprintf("median %g genes per chromosome after the expression cutoff (%d genes, %d chromosomes); smoothing window %d genes",
                   med, length(chr), length(n), window_length)
    flog.info(sprintf("HMM gene density check: %s", msg))
    short <- names(n)[n < window_length / 2]
    if (length(short) && med >= window_length) {
        flog.warn(sprintf("chromosomes with fewer than %g genes (HMM calls there are less reliable): %s", window_length / 2,
                          paste(sprintf("%s (%d)", short, n[short]), collapse=", ")))
    }
    if (med < window_length) {
        txt <- sprintf(paste0("Too few genes for the HMM: %s. Most chromosomes are shorter than one smoothing window, so the ",
                              "HMM would call noise as CNVs. Use more reads per cell (e.g. larger spatial bins, or pooled cells), ",
                              "or run with HMM=FALSE; set HMM_gene_density_check='warn' to proceed anyway."), msg)
        if (action == "stop") {
            flog.error(txt)
            stop(txt)
        }
        flog.warn(txt)
    } else if (med < 2 * window_length) {
        flog.warn(sprintf("Low gene density for the HMM: %s. CNV calls may be noisy; compare them with the expression heatmap.", msg))
    }
    invisible(list(median=med, per_chr=n))
}


# Empirical false-call rate of the i6 HMM. Reference cells carry no CNVs, so the HMM, run with the emission model it
# will use on the observations, on random groups of reference cells (drawn within each reference group, sized like the
# observation subclusters) should call (almost) everything neutral; the fraction of genes x groups called non-neutral
# estimates how often the run's calls are noise. Random groups are an optimistic null (they lack the shared structure of
# real subclusters), so the thresholds are set low. Stop (or warn) above thresholds[2]; warn above thresholds[1].
# The random number state is restored afterwards so the rest of the run is unchanged.
.check_hmm_null_fpr <- function(infercnv_obj, t, action=c("stop", "warn"), thresholds=c(0.05, 0.15), max_groups=60,
                                min_groups=5, out_dir=NULL,
                                cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)) {
    action <- match.arg(action)
    ref_groups <- infercnv_obj@reference_grouped_cell_indices
    if (length(ref_groups) == 0) {
        flog.info("HMM null false-call check: no reference cells, skipped")
        return(invisible(NULL))
    }
    ref_cells <- unlist(ref_groups)
    if (!is.null(infercnv_obj@tumor_subclusters)) {
        sub <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
        sizes <- lengths(sub[!vapply(sub, function(ix) all(ix %in% ref_cells), logical(1))])
    } else {
        sizes <- lengths(infercnv_obj@observation_grouped_cell_indices)
    }
    if (length(sizes) == 0) {
        return(invisible(NULL))
    }
    has_seed <- exists(".Random.seed", envir=globalenv())
    if (has_seed) {
        old_seed <- get(".Random.seed", envir=globalenv())
        on.exit(assign(".Random.seed", old_seed, envir=globalenv()))
    }
    draw <- function(sizes) {
        groups <- list()
        for (g in names(ref_groups)) {
            pool <- ref_groups[[g]]
            pool <- pool[sample.int(length(pool))]
            while (length(pool) >= min(sizes) && length(groups) < max_groups) {
                k <- min(sizes[sample.int(length(sizes), 1)], length(pool))
                groups[[length(groups) + 1]] <- pool[seq_len(k)]
                pool <- pool[-seq_len(k)]
            }
        }
        groups
    }
    set.seed(1)
    groups <- draw(sizes)
    sized <- "sized like the observation subclusters"
    if (length(groups) < min_groups) {
        ## observation groups larger than the reference pool (e.g. analysis_mode='samples'): smaller groups, with the
        ## emission sd for their size
        cap <- max(floor(sum(lengths(ref_groups)) / (2 * min_groups)), 1)
        flog.info(sprintf("HMM null false-call check: observation groups larger than the reference pool; using groups of up to %d reference cells", cap))
        groups <- draw(pmin(sizes, cap))
        sized <- "smaller than the observation groups"
    }
    if (length(groups) < min_groups) {
        flog.info(sprintf("HMM null false-call check: only %d groups of reference cells the size of the observation subclusters, skipped",
                          length(groups)))
        return(invisible(NULL))
    }
    HMM_info <- .get_HMM(cnv_mean_sd, t)
    gene_order <- infercnv_obj@gene_order
    chr_idx <- split(seq_len(nrow(gene_order)), factor(gene_order$chr, levels=unique(gene_order$chr)))
    called <- vapply(groups, function(cells) {
        ep <- .get_state_emission_params(length(cells), cnv_mean_sd, cnv_level_to_mean_sd_fit)
        st <- unlist(lapply(chr_idx, function(gi) {
            v <- rowMeans(infercnv_obj@expr.data[gi, cells, drop=FALSE])
            Viterbi.dthmm.adj(HiddenMarkov::dthmm(v, HMM_info[['state_transitions']], HMM_info[['delta']], "norm", ep))
        }))
        mean(st != 3)
    }, numeric(1))
    fpr <- mean(called)
    msg <- sprintf("%.1f%% of genes called as CNVs in %d random groups of reference cells (median %d cells, %s; worst group %.1f%%)",
                   100 * fpr, length(groups), as.integer(stats::median(lengths(groups))), sized, 100 * max(called))
    flog.info(sprintf("HMM null false-call check: %s", msg))
    if (!is.null(out_dir)) {
        utils::write.table(data.frame(groups=length(groups), median_group_cells=stats::median(lengths(groups)),
                                      null_false_call_rate=fpr, worst_group=max(called),
                                      warn_above=thresholds[1], stop_above=thresholds[2]),
                           file.path(out_dir, "HMM_null_false_call_rate.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
    if (fpr > thresholds[2]) {
        txt <- sprintf(paste0("HMM calls would be unreliable: %s, above %g%%. The data are too noisy for the HMM to separate ",
                              "copy-number states. Use more reads per cell (e.g. larger spatial bins, or pooled cells), or run ",
                              "with HMM=FALSE; set HMM_null_fpr_check='warn' to proceed anyway."), msg, 100 * thresholds[2])
        if (action == "stop") {
            flog.error(txt)
            stop(txt)
        }
        flog.warn(txt)
    } else if (fpr > thresholds[1]) {
        flog.warn(sprintf("Elevated HMM false-call rate: %s. Treat small or low-amplitude calls with caution.", msg))
    }
    invisible(list(rate=fpr, per_group=called))
}


#' @title .get_HMM
#'
#' @description retrieves parameters for the i6 HMM including state transition and state emission probabilities.
#'
#' @param cnv_mean_sd list containing the N(mean,sd) values per cnv state level
#'
#' @param t  the probability for transitioning to a different state.
#'
#' @return HMM_info list
#'
#' @noRd


.get_HMM <- function(cnv_mean_sd, t) {

    ##                      states: 0,   0.5,    1,     1.5,     2,    3
    state_transitions = matrix( c(1-5*t,  t,     t,      t,      t,    t,
                                  t,     1-5*t,  t,      t,      t,    t,
                                  t,      t,     1-5*t,  t,      t,    t,
                                  t,      t,     t,    1-5*t,    t,    t,
                                  t,      t,     t,      t,    1-5*t,  t,
                                  t,      t,     t,      t,      t,  1-5*t),
                               byrow=TRUE,
                               nrow=6)

    delta=c(t,      t,     1-5*t,    t,    t,   t) # more likely normal,
    
    state_emission_params = list(mean=c(cnv_mean_sd[["cnv:0.01"]]$mean,
                                        cnv_mean_sd[["cnv:0.5"]]$mean,
                                        cnv_mean_sd[["cnv:1"]]$mean,
                                        cnv_mean_sd[["cnv:1.5"]]$mean,
                                        cnv_mean_sd[["cnv:2"]]$mean,
                                        cnv_mean_sd[["cnv:3"]]$mean),

                                 sd=c(cnv_mean_sd[["cnv:0.01"]]$sd,
                                      cnv_mean_sd[["cnv:0.5"]]$sd,
                                      cnv_mean_sd[["cnv:1"]]$sd,
                                      cnv_mean_sd[["cnv:1.5"]]$sd,
                                      cnv_mean_sd[["cnv:2"]]$sd,
                                      cnv_mean_sd[["cnv:3"]]$sd) )
    
    
    
    HMM_info = list(state_transitions=state_transitions,
                    delta=delta,
                    state_emission_params=state_emission_params)
    
    return(HMM_info)
}


#' @title predict_CNV_via_HMM_on_indiv_cells
#'
#' @description predict CNV levels at the individual cell level, using the i6 HMM
#'
#' @param infercnv_obj infercnv object
#'
#' @param cnv_mean_sd (optional, by default automatically computed based in the infercnv_obj@.hspike object)
#'
#' @param t HMM alt state transition probability (default=1e-6)
#'
#' @return infercnv_obj where the infercnv_obj@expr.data are replaced with the HMM state assignments.
#'
#' @keywords internal
#' @noRd
#'

predict_CNV_via_HMM_on_indiv_cells  <- function(infercnv_obj, cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike), t=1e-6) {
    
    flog.info("predict_CNV_via_HMM_on_indiv_cells()")
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

        
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(seq_len(ncol(expr.data)), function(cell_idx) {
            
            gene_expr_vals = as.vector(expr.data[chr_gene_idx,cell_idx])
            
            hmm <- HiddenMarkov::dthmm(x=gene_expr_vals,
                                       Pi=HMM_info[['state_transitions']],
                                       delta=HMM_info[['delta']],
                                       distn="norm",
                                       pm=HMM_info[['state_emission_params']])
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}


#' @title predict_CNV_via_HMM_on_tumor_subclusters
#'
#' @description predict CNV levels at the tumor subcluster level, using the i6 HMM
#'
#' @param infercnv_obj infercnv object
#'
#' @param cnv_mean_sd (optional, by default automatically computed based in the infercnv_obj@.hspike object)
#'
#' @param cnv_level_to_mean_sd_fit (optional, by default automatically computed based on get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)
#' 
#' @param t HMM alt state transition probability (default=1e-6)
#'
#' @return infercnv_obj where the infercnv_obj@expr.data are replaced with the HMM state assignments.
#'
#' @keywords internal
#' @noRd
#'

predict_CNV_via_HMM_on_tumor_subclusters  <- function(infercnv_obj,
                                                      cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                      cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                      t=1e-6
                                                      ) {
    
    


    flog.info("predict_CNV_via_HMM_on_tumor_subclusters")

    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("No subclusters defined, so instead running on whole samples")
        # arguments by name: the second positional argument is cluster_by_groups
        return(predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, cluster_by_groups=TRUE, cnv_mean_sd=cnv_mean_sd,
                                                          cnv_level_to_mean_sd_fit=cnv_level_to_mean_sd_fit, t=t));
    }
    
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    tumor_subclusters <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
    
    ## add the normals, so they get predictions too:
    #tumor_subclusters <- c(tumor_subclusters, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_subclusters, function(tumor_subcluster_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx,drop=FALSE])
            
            num_cells = length(tumor_subcluster_cells_idx)

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)

            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)

            hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}



predict_CNV_via_HMM_on_tumor_subclusters_per_chr  <- function(infercnv_obj,
                                                              subclusters_per_chr,
                                                              cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                              cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                              t=1e-6
                                                              ) {
    

    flog.info("predict_CNV_via_HMM_on_tumor_subclusters_per_chr")

    if (is.null(subclusters_per_chr)) {
        flog.warn("No subclusters defined, so instead running on whole samples")
        # arguments by name: the second positional argument is cluster_by_groups
        return(predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, cluster_by_groups=TRUE, cnv_mean_sd=cnv_mean_sd,
                                                          cnv_level_to_mean_sd_fit=cnv_level_to_mean_sd_fit, t=t));
    }
    
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    tumor_subclusters <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
    
    ## add the normals, so they get predictions too:
    #tumor_subclusters <- c(tumor_subclusters, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(subclusters_per_chr[[chr]], function(subcluster_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,subcluster_cells_idx,drop=FALSE])
            
            num_cells = length(subcluster_cells_idx)

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
            
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    flog.info("-done predicting CNV based on per chromosome subclusters")

    flog.info("-calculating initial tumor subclusters CNV consensus based on per chromosome predictions")

    ret = get_predicted_CNV_regions(infercnv_obj, by="subcluster", mode="expression")
    for (i in seq_along(ret)) {
        for (region in ret[[i]]$gene_regions) {
            if (length(unique(region$state)) != 1) {
                flog.error("error, all states are not equal")
            }
            else {
                infercnv_obj@expr.data[region$gene, ret[[i]]$cells] = rep(region$state[1], length(region$gene)*length(ret[[i]]$cells))
            }
        }    
    }

    return(infercnv_obj)
    
}


#' @title predict_CNV_via_HMM_on_whole_tumor_samples
#'
#' @description predict CNV levels at the tumor sample level, using the i6 HMM
#'
#' @param infercnv_obj infercnv object
#'
#' @param cnv_mean_sd (optional, by default automatically computed based in the infercnv_obj@.hspike object)
#'
#' @param cnv_level_to_mean_sd_fit (optional, by default automatically computed based on get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)
#'
#' @param t HMM alt state transition probability (default=1e-6)
#'
#' @return infercnv_obj where the infercnv_obj@expr.data are replaced with the HMM state assignments.
#'
#' @keywords internal
#' @noRd
#'


predict_CNV_via_HMM_on_whole_tumor_samples  <- function(infercnv_obj,
                                                        cluster_by_groups,
                                                        cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                                        cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike),
                                                        t=1e-6
                                                        ) {
    
    
    flog.info("predict_CNV_via_HMM_on_whole_tumor_samples")
    
    HMM_info  <- .get_HMM(cnv_mean_sd, t)
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    ## add the normals, so they get predictions too:
    if (cluster_by_groups == TRUE) {
        tumor_samples = c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    } else {
        tumor_samples = c(all_observations = unlist(infercnv_obj@observation_grouped_cell_indices), infercnv_obj@reference_grouped_cell_indices)
    }
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_samples, function(tumor_sample_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_sample_cells_idx,drop=FALSE])
            
            num_cells = length(tumor_sample_cells_idx)

            state_emission_params <- .get_state_emission_params(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit)
                        
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       state_emission_params)
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_sample_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}


#' @title .get_state_emission_params
#'
#' @description Given a specified number of cells, determines the standard deviation for each of the cnv states
#' based on the linear model fit.
#'
#' @param num_cells  number of cells in the tumor subcluster
#'
#' @param cnv_mean_sd list of cnv mean,sd values
#'
#' @param cnv_level_to_mean_sd_fit  linear model that was fit for each cnv state level
#'
#' @param plot boolean (default=FALSE)
#'
#' @noRd


.get_state_emission_params <- function(num_cells, cnv_mean_sd, cnv_level_to_mean_sd_fit, plot=FALSE) {
        
    for (cnv_level in names(cnv_mean_sd)) {
        sd <- .hspike_sd_at(cnv_level_to_mean_sd_fit[[cnv_level]], num_cells)
        cnv_mean_sd[[cnv_level]]$sd <- sd
    }
    
    if (plot) {
        .plot_cnv_mean_sd_for_num_cells(num_cells, cnv_mean_sd)
    }

    state_emission_params = list(mean=c(cnv_mean_sd[["cnv:0.01"]]$mean,
                                        cnv_mean_sd[["cnv:0.5"]]$mean,
                                        cnv_mean_sd[["cnv:1"]]$mean,
                                        cnv_mean_sd[["cnv:1.5"]]$mean,
                                        cnv_mean_sd[["cnv:2"]]$mean,
                                        cnv_mean_sd[["cnv:3"]]$mean),
                                 
                                 sd=c(cnv_mean_sd[["cnv:0.01"]]$sd,
                                      cnv_mean_sd[["cnv:0.5"]]$sd,
                                      cnv_mean_sd[["cnv:1"]]$sd,
                                      cnv_mean_sd[["cnv:1.5"]]$sd,
                                      cnv_mean_sd[["cnv:2"]]$sd,
                                      cnv_mean_sd[["cnv:3"]]$sd) )
    
    
    return(state_emission_params)
}



#' @title .plot_cnv_mean_sd_for_num_cells
#'
#' @description helper function for plotting the N(mean,sd) for each of the cnv states
#' The plot is written to a file 'state_emissions.{num_cells}.pdf
#' 
#' @param num_cells number of cells. Only used to encode into the filename.
#'
#' @param cnv_mean_sd list containing the N(mean,sd) for each of the cnv states
#'
#' @return None
#'
#' @noRd

.plot_cnv_mean_sd_for_num_cells <- function(num_cells, cnv_mean_sd) {

    pdf(sprintf("state_emissions.%d-cells.pdf", num_cells))
    
    # p = ggplot(data.frame(x=c(0,3)), aes(x)) +
    p = ggplot(data.frame(x=c(0,3)), aes_string(x='x')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:0.01"]]$mean,'sd'=cnv_mean_sd[["cnv:0.01"]]$sd),
                      aes(colour='cnv:0')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:0.5"]]$mean,'sd'=cnv_mean_sd[["cnv:0.5"]]$sd),
                      aes(colour='cnv:0.5')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:1"]]$mean,'sd'=cnv_mean_sd[["cnv:1"]]$sd),
                      aes(colour='cnv:1')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:1.5"]]$mean,'sd'=cnv_mean_sd[["cnv:1.5"]]$sd),
                      aes(colour='cnv:1.5')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:2"]]$mean,'sd'=cnv_mean_sd[["cnv:2"]]$sd),
                      aes(colour='cnv:2')) +
        stat_function(fun=dnorm,
                      args=list('mean'=cnv_mean_sd[["cnv:3"]]$mean,'sd'=cnv_mean_sd[["cnv:3"]]$sd),
                      aes(colour='cnv:3')) +
        scale_colour_manual("Function",
                            values=c("blue", "red", "green", "magenta", "orange", "cyan"),
                            breaks=c("cnv:0","cnv:0.5", "cnv:1", "cnv:1.5", "cnv:2", "cnv:3"))
    
    plot(p)

    dev.off()

}



#' @keywords internal
#' @noRd
#' 


.compare_obs_vs_fit_sc_cnv_sd <- function(infercnv_obj,
                                          cnv_mean_sd=get_spike_dists(infercnv_obj@.hspike),
                                          cnv_level_to_mean_sd_fit=get_hspike_cnv_mean_sd_trend_by_num_cells_fit(infercnv_obj@.hspike)) {
    
    cnv_levels = names(cnv_mean_sd)

    df = do.call(rbind, lapply(cnv_levels, function(cnv_level) {
        sd_obs=cnv_mean_sd[[cnv_level]]$sd
        sd_pred=.hspike_sd_at(cnv_level_to_mean_sd_fit[[cnv_level]], 1)
        return(data.frame(cnv=cnv_level, sd_obs=sd_obs, sd_pred=sd_pred))

    }))
    
    return(df)
}



#' @title get_predicted_CNV_regions
#'
#' @description Given the infercnv_obj containing the HMM state assignments in the expr.data slot,
#' retrieves a list of CNV regions.
#'
#' @param infercnv_obj infercnv object
#'
#' @param by options("consensus", "subcluster", "cell"), determines the granularity at which to report
#'           the CNV regions.  Ideally, set to the same level at which the HMM predictions were performed.
#'
#' @return cnv_regions list
#'
#' @keywords internal
#' @noRd
#'

get_predicted_CNV_regions <- function(infercnv_obj,
                                      by=c("consensus", "subcluster", "cell"),
                                      mode=c("expression", "allele")) {
    by = match.arg(by)
    mode = match.arg(mode)

    if (mode == "expression") {
        expr.data = infercnv_obj@expr.data
        gene_order = infercnv_obj@gene_order
    } else { # "allele"
        expr.data = infercnv_obj@allele.expr.data
        gene_order = infercnv_obj@allele.gene_order
    }

    flog.info(sprintf("get_predicted_CNV_regions(%s)", by)) 
    
    cell_groups = NULL

    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("get_predicted_CNV_regions() - no subclusters defined, resetting reporting mode to consensus")
        by <- "consensus"
    }
    
    if (by == "consensus") {
        # cell_groups = infercnv_obj@observation_grouped_cell_indices
        cell_groups = c(infercnv_obj@reference_grouped_cell_indices, infercnv_obj@observation_grouped_cell_indices)
    } else if (by == "subcluster") {
        cell_groups = unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
    } else if (by == "cell") {
        # cell_groups = lapply(unlist(infercnv_obj@observation_grouped_cell_indices), function(x) x) 
        cell_groups = c(unlist(infercnv_obj@reference_grouped_cell_indices, use.names=FALSE), unlist(infercnv_obj@observation_grouped_cell_indices, use.names=FALSE))
        names(cell_groups) = colnames(expr.data)[cell_groups]
    }
    else {
        stop("Error, shouldn't get here ... bug")
    }
    
    cnv_regions = list()

    cnv_counter_start = 0
    for (cell_group_name in names(cell_groups)) {
        
        cell_group = cell_groups[[cell_group_name]]
        #flog.info(sprintf("cell group %s -> %s", cell_group_name, cell_group))

        flog.info(sprintf("-processing cell_group_name: %s, size: %d", cell_group_name, length(cell_group)))
                
        cell_group_mtx = expr.data[,cell_group,drop=FALSE]
        cell_group_names = colnames(cell_group_mtx)

        state_consensus <- .get_state_consensus(cell_group_mtx)
        
        names(state_consensus) <- rownames(cell_group_mtx)
        cnv_gene_regions <- .define_cnv_gene_regions(state_consensus, gene_order, cnv_counter_start)
        cnv_ranges <- .get_cnv_gene_region_bounds(cnv_gene_regions)

        consensus_state_list = list(cell_group_name=cell_group_name,
                                    cells=cell_group_names,
                                    gene_regions=cnv_gene_regions,
                                    cnv_ranges=cnv_ranges)
        
        cnv_regions[[length(cnv_regions)+1]] = consensus_state_list

        cnv_counter_start = cnv_counter_start + length(cnv_gene_regions)
        
    }
    
    return(cnv_regions)
    
}


#' @title generate_cnv_region_reports
#'
#' @description writes the CNV region report files
#'
#' @param infercnv_obj infercnv object
#'
#' @param output_filename_prefix  prefix for output filename
#' 
#' @param out_dir output directory for report files to be written
#'
#' @param ignore_neutral_state  numeric value representing the neutral state, which should be excluded from reporting (default: NA)
#' 
#' @param by options("consensus", "subcluster", "cell"), determines the granularity at which to report
#'           the CNV regions.  Ideally, set to the same level at which the HMM predictions were performed.
#'
#' @param mode options("expression", "allele"), determines if to generate the reports on the expression based results or allele based results.
#'
#' @return None
#'
#' @keywords internal
#' @noRd
#'



generate_cnv_region_reports <- function(infercnv_obj,
                                        output_filename_prefix,
                                        out_dir,
                                        ignore_neutral_state=NA,
                                        by=c("consensus", "subcluster", "cell"),
                                        mode=c("expression", "allele")) {

    mode = match.arg(mode)

    cnv_regions <- get_predicted_CNV_regions(infercnv_obj, by, mode)

    ## cell clusters defined.
    cell_clusters_outfile = paste(out_dir, paste0(output_filename_prefix, ".cell_groupings"), sep="/")

    cell_clusters_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        cells = x$cells
        ret_df = data.frame(cell_group_name=cell_group_name, cell=cells)
        return(ret_df)
    })
    
    cell_clusters_df = do.call(rbind, cell_clusters_df)
    flog.info(sprintf("-writing cell clusters file: %s", cell_clusters_outfile))
    write.table(cell_clusters_df, file=cell_clusters_outfile, row.names=FALSE, quote=FALSE, sep="\t")
    
    ## regions DF:
    regions_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_regions.dat"), sep="/")
    
    regions_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        cnv_ranges = x$cnv_ranges
        ret_df = cbind(cell_group_name=cell_group_name, cnv_ranges)
        return(ret_df)
    })

    regions_df = do.call(rbind, regions_df)

    ## remove the neutral calls:
    if (! is.na(ignore_neutral_state)) {
        regions_df = regions_df[regions_df$state != ignore_neutral_state, ]
    }
    flog.info(sprintf("-writing cnv regions file: %s", regions_outfile)) 
    write.table(regions_df, regions_outfile, row.names=FALSE, sep="\t", quote=FALSE)


    ## write the per-gene reports of cnv:
    if (by == "cell") {
        flog.warn("Note, HMM reporting is being done by 'cell', so this may use more memory, write more info to disk, take more time, ...")
    }
    gene_cnv_df = lapply(cnv_regions, function(x) {
        cell_group_name = x$cell_group_name
        gene_region_list = x$gene_regions
        gene_region_names = names(gene_region_list)
        gene_region_df = lapply(gene_region_names, function(gene_region_name) {
            df = gene_region_list[[gene_region_name]]
            df = cbind(cell_group_name, gene_region_name, df)
            rownames(df) <- NULL
            return(df)
        })
        gene_region_df = do.call(rbind, gene_region_df)

        return(gene_region_df)
    })
    gene_cnv_df = do.call(rbind, gene_cnv_df)
    if (! is.na(ignore_neutral_state)) {
        gene_cnv_df = gene_cnv_df[gene_cnv_df$state != ignore_neutral_state, ]
    }
    
    ## write output file:
    gene_cnv_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_genes.dat"), sep="/") 
    flog.info(sprintf("-writing per-gene cnv report: %s", gene_cnv_outfile))
    write.table(gene_cnv_df, gene_cnv_outfile, row.names=FALSE, sep="\t", quote=FALSE)

    ## write file containing all genes that were leveraged in the predictions:
    gene_order_outfile = paste(out_dir, paste0(output_filename_prefix, ".genes_used.dat"), sep="/")
    flog.info(sprintf("-writing gene ordering info: %s", gene_order_outfile))
    if (mode == "expression") {
        write.table(infercnv_obj@gene_order, file=gene_order_outfile, quote=FALSE, sep="\t")
    } else { # "allele"
        write.table(infercnv_obj@allele.gene_order, file=gene_order_outfile, quote=FALSE, sep="\t")
    }
    
    
    return
    
}



#' @title adjust_genes_regions_report
#'
#' @description Updates the CNV region report files
#'
#' @param infercnv_obj infercnv object
#'
#' @param input_filename_prefix  prefix for input filename
#' 
#' @param output_filename_prefix  prefix for output filename
#' 
#' @param out_dir output directory for report files to be written
#' 
#' @return None
#'
#' @keywords internal
#' @noRd
#'

adjust_genes_regions_report <- function(hmm.infercnv_obj,
                                        input_filename_prefix,
                                        output_filename_prefix,
                                        out_dir) {
    #####################
    # Load Data
    #####################
    # Get the paths to the files 
    pred_cnv_genes_path <- paste(out_dir, paste0(input_filename_prefix, ".pred_cnv_genes.dat"), sep="/")
    pred_cnv_regions_path <- paste(out_dir, paste0(input_filename_prefix, ".pred_cnv_regions.dat"), sep="/")
    
    # Check if files exist 
    if (file.exists(pred_cnv_genes_path)){
        # Read in the files 
        pred_cnv_genes <- read.table(file = pred_cnv_genes_path, header = T, sep = "\t")
    } else {
        # If the file cant be found, throw the error message 
        error_message <- paste("Cannot find and adjust the following file.", pred_cnv_genes_path)
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    
    # Check if files exist 
    if (file.exists(pred_cnv_genes_path)){
        # Read in the files 
        pred_cnv_regions <- read.table(file = pred_cnv_regions_path, header = T, sep = "\t")
    } else{
        # If the file cant be found, throw the error message 
        error_message <- paste("Cannot find and adjust the following file.", pred_cnv_regions_path)
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    
    #####################
    # Process the data 
    #####################
    # CNV regions that are in the the new data set after removal of cnv's 
    new_regions <- sapply(hmm.infercnv_obj@cell_gene, function(i) i$cnv_regions)
    
    # Subset the file to only include the cnv that were not removed 
    #-----------------
    ## Genes file
    #-----------------
    ids <- which(pred_cnv_genes$gene_region_name %in% new_regions)
    new_pred_cnv_genes <- pred_cnv_genes[ids, ]
    
    # Get the new states predicted by the bayesian model and add them to the data frame 
    new_states <- sapply(new_pred_cnv_genes$gene_region_name, function(i) {
        hmm.infercnv_obj@cell_gene[[which(new_regions %in% i)]]$State
    })
    # Assign the new states to the state column in the data frame 
    new_pred_cnv_genes$state <- new_states
    
    # write the new file 
    gene_cnv_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_genes.dat"), sep="/") 
    write.table(new_pred_cnv_genes, gene_cnv_outfile, row.names=FALSE, sep="\t", quote=FALSE)
    
    #-----------------
    ## Regions file
    #-----------------
    ids <- which(pred_cnv_regions$cnv_name %in% new_regions)
    new_pred_cnv_regions <- pred_cnv_regions[ids, ]
    
    # Get the new states predicted by the bayesian model and add them to the data frame 
    new_states <- sapply(new_pred_cnv_regions$cnv_name, function(i) {
        hmm.infercnv_obj@cell_gene[[which(new_regions %in% i)]]$State
    })
    # Assign the new states to the state column in the data frame 
    new_pred_cnv_regions$state <- new_states
    # write the new file 
    regions_cnv_outfile = paste(out_dir, paste0(output_filename_prefix, ".pred_cnv_regions.dat"), sep="/") 
    write.table(new_pred_cnv_regions, regions_cnv_outfile, row.names=FALSE, sep="\t", quote=FALSE)
}



#' @title .get_state_consensus
#'
#' @description gets the state consensus for each gene in the input matrix
#'
#' @param cell_group_matrix  matrix of [genes,cells] for which to apply consensus operation at gene level across cells.
#'
#' @return vector containing consensus state assignments for genes.
#'
#' @noRd

.get_state_consensus <- function(cell_group_matrix) {

    consensus  = apply(cell_group_matrix, 1, function(x) {
        t = table(x)
        names(t)[order(t, decreasing=TRUE)[1]]
    })

    consensus <- as.numeric(consensus)
    
    return(consensus)
}


#' @title .define_cnv_gene_regions
#'
#' @description Given the state consensus vector and gene order info, defines cnv regions
#' based on consistent ordering and cnv state 
#'
#' @param state_consensus state consensus vector
#'
#' @param gene_order the infercnv_obj@gene_order info
#'
#' @param cnv_region_counter number x where counting starts at x+1, used to provide unique region names.
#'
#' @return regions  list containing the cnv regions defined.
#'
#' @noRd

.define_cnv_gene_regions <- function(state_consensus, gene_order, cnv_region_counter) {

    regions = list()

    gene_names = rownames(gene_order)
    
    chrs = unique(gene_order$chr)
    for (chr in chrs) {
        gene_idx = which(gene_order$chr==chr)
        if (length(gene_idx) < 2) { next }
        
        chr_states = state_consensus[gene_idx]
        prev_state = chr_states[1]
        ## pos_begin = paste(gene_order[gene_idx[1],,drop=TRUE], collapse=",")
        pos_begin = gene_order[gene_idx[1],,drop=TRUE]
        
        cnv_region_counter = cnv_region_counter + 1

        cnv_region_name = sprintf("%s-region_%d", chr, cnv_region_counter)
        current_cnv_region = data.frame(state=prev_state,
                                        gene=gene_names[gene_idx[1]],
                                        chr=pos_begin$chr,
                                        start=pos_begin$start,
                                        end=pos_begin$stop) 
        regions[[cnv_region_name]] = current_cnv_region
        
        for (i in seq(2,length(gene_idx))) {
            state = chr_states[i]
            pos_end = gene_order[gene_idx[i],,drop=TRUE]
            next_gene_entry = data.frame(state=state,
                                         gene=gene_names[gene_idx[i]],
                                         chr=pos_end$chr,
                                         start=pos_end$start,
                                         end=pos_end$stop)
            
            if (state != prev_state) {
                ## state transition
                ## start new cnv region
                cnv_region_counter = cnv_region_counter + 1
                cnv_region_name = sprintf("%s-region_%d", chr, cnv_region_counter)
                regions[[cnv_region_name]] = next_gene_entry
            } else {
                ## append gene to current cnv region
                regions[[cnv_region_name]] = rbind(regions[[cnv_region_name]], next_gene_entry)
            }

            prev_state = state
        }
        
    }

    return(regions)
}


#' @title .get_cnv_gene_region_bounds
#'
#' @description  Given the cnv regions list, defines a data table containing
#'               the cnv region name, state, chr, start, and end value.
#'
#' @param cnv_gene_regions
#'
#' @return data.frame containing the cnv region summary table
#'
#' @noRd

.get_cnv_gene_region_bounds <- function(cnv_gene_regions) {

    bounds = do.call(rbind, lapply(names(cnv_gene_regions), function(x) {
        cnv_name = x
        df = cnv_gene_regions[[cnv_name]]
        state = df$state[1]
        chr = df$chr[1]
        start = min(df$start)
        end = max(df$end)
        
        return(data.frame(cnv_name, state, chr, start, end))
    }) )

    rownames(bounds) <- NULL
    
    return(bounds)
}


#' @title Viterbi.dthmm.adj
#'
#' @description Viterbi method extracted from the HiddenMarkov package and modified
#'              to use our scoring system.
#'
#' @param HiddenMarkov object
#'
#' @return vector containing the viterbi state assignments
#'
#' @noRd

Viterbi.dthmm.adj <- function (object, ...){
    x <- object$x

    if (length(x) < 2) {
        ## not enough run a trace on
        return(which.max(object$delta)); # neutral state (3 for i6, 2 for i3)
    }
    
    dfunc <- HiddenMarkov:::makedensity(object$distn)
    n <- length(x)
    m <- nrow(object$Pi) # transition matrix
    nu <- matrix(NA, nrow = n, ncol = m)  # scoring matrix
    y <- rep(NA, n) # final trace
    pseudocount = 1e-20
    
    emissions <- matrix(NA, nrow = n, ncol = m) 
    
    ## ###############################
    ## restrict to constant variance to avoid nonsensical results:

    ## object$pm$sd = max(object$pm$sd)
    object$pm$sd = median(object$pm$sd) # max is too high
    
    ## ###############################    


    ## init first row

    emission <- pnorm(abs(x[1]-object$pm$mean)/object$pm$sd, log.p=TRUE, lower.tail=FALSE)
    emission <- 1 / (-1 * emission)
    emission <- emission / sum(emission)
    
    emissions[1,] <- log(emission)
    
    nu[1, ] <- log(object$delta) + # start probabilities
        emissions[1,]
    
    
    #nu[1, ] <- log(object$delta) + # start probabilities
    #    dfunc(x=x[1], # mean expr val for gene_1
    #          object$pm, # parameters (mean, sd) for norm dist
    #          HiddenMarkov:::getj(object$pn, 1), # NULL value
    #          log=TRUE) # returns p-values as log(p)
    
    logPi <- log(object$Pi) # convert transition matrix to log(p)
    
    for (i in 2:n) {
        
        matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
        
        #nu[i, ] <- apply(matrixnu + logPi, 2, max) +
        #              dfunc(x=x[i], object$pm, getj(object$pn, i),
        #                    log=TRUE)

        
        emission <- pnorm(abs(x[i]-object$pm$mean)/object$pm$sd, log.p=TRUE, lower.tail=FALSE)
        emission <- 1 / (-1 * emission)
        emission <- emission / sum(emission)
        
        emissions[i, ] <- log(emission)
        
        nu[i, ] <- apply(matrixnu + logPi, 2, max) + emissions[i, ] 
                
    }
    if (any(nu[n, ] == -Inf)) 
        stop("Problems With Underflow")


    ## traceback
    y[n] <- which.max(nu[n, ])

    for (i in seq(n - 1, 1, -1))
        y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])

    return(y)
}


#' @title assign_HMM_states_to_proxy_expr_vals
#'
#' @description Replaces the HMM state assignments with the cnv levels they represent.
#' 
#' @param infercnv_obj infercnv object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

assign_HMM_states_to_proxy_expr_vals <- function(infercnv_obj) {

    expr.data = infercnv_obj@expr.data

    expr.data[expr.data == 1] <- 0
    expr.data[expr.data == 2] <- 0.5
    expr.data[expr.data == 3] <- 1
    expr.data[expr.data == 4] <- 1.5
    expr.data[expr.data == 5] <- 2
    expr.data[expr.data == 6] <- 3

    infercnv_obj@expr.data <- expr.data

    return(infercnv_obj)

}

