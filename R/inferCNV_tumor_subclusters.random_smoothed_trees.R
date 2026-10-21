

define_signif_tumor_subclusters_via_random_smooothed_trees <- function(infercnv_obj, p_val, hclust_method, cluster_by_groups, window_size=101,
                                                                       max_recursion_depth=3, min_cluster_size_recurse=10) {
    
    ## the state of the infercnv object here should be:
    ## log transformed
    ## but *NOT* smoothed.
    ## TODO: -include check for smoothed property so will not run this if already smoothed.
        
    infercnv_copy = infercnv_obj  ## don't want to change the original data .... just want to add subcluster info.
    
    flog.info(sprintf("define_signif_tumor_subclusters(p_val=%g", p_val))

    infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=TRUE)  # important, remove normal from tumor before testing clusters.

    ## must treat normals same way!
    tumor_groups = list()
    if (cluster_by_groups) {
        tumor_groups <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    }
    else {
        # if(length(infercnv_obj@reference_grouped_cell_indices) > 0) {
        tumor_groups <- c(list(all_observations=unlist(infercnv_obj@observation_grouped_cell_indices, use.names=FALSE)), infercnv_obj@reference_grouped_cell_indices)
        # }
        # else {
        #     tumor_groups <- list(all_observations=unlist(infercnv_obj@observation_grouped_cell_indices, use.names=FALSE))
        # }
    }

    res = list()
    
    for (tumor_group in names(tumor_groups)) {
        
        flog.info(sprintf("define_signif_tumor_subclusters(), tumor: %s", tumor_group))
        
        tumor_group_idx <- tumor_groups[[ tumor_group ]]
        names(tumor_group_idx) = colnames(infercnv_obj@expr.data)[tumor_group_idx]
        tumor_expr_data <- infercnv_obj@expr.data[,tumor_group_idx]
                
        tumor_subcluster_info <- .single_tumor_subclustering_smoothed_tree(tumor_group, tumor_group_idx, tumor_expr_data, p_val, hclust_method, window_size,
                                                                           max_recursion_depth, min_cluster_size_recurse)
        
        res$hc[[tumor_group]] <- tumor_subcluster_info$hc
        res$subclusters[[tumor_group]] <- tumor_subcluster_info$subclusters
        
    }
    
    infercnv_copy@tumor_subclusters <- res
    
    if (! is.null(infercnv_copy@.hspike)) {
        flog.info("-mirroring for hspike")

        # as for the other partition methods (see define_signif_tumor_subclusters), the hspike is not subclustered,
        # which would lower the expected noise level of the states. (This call previously passed its arguments by
        # position, shifted by one, so the hspike was subclustered with window_size = 3.)
        infercnv_copy@.hspike <- define_signif_tumor_subclusters(infercnv_copy@.hspike,
                                                                 cluster_by_groups = TRUE,
                                                                 partition_method = "none")[[1]]
    }
    
    
    return(infercnv_copy)
}



.single_tumor_subclustering_smoothed_tree <- function(tumor_name, tumor_group_idx, tumor_expr_data, p_val, hclust_method, window_size,
                                                      max_recursion_depth, min_cluster_size_recurse) {


    tumor_subcluster_info = list()

    ## smooth and median-center
    sm_tumor_expr_data = apply(tumor_expr_data, 2, caTools::runmean, k=window_size)
    #sm_tumor_expr_data = scale(sm_tumor_expr_data, center=TRUE, scale=FALSE)
    sm_tumor_expr_data = .center_columns(sm_tumor_expr_data, 'median')
        
    
    hc <- hclust(parallelDist(t(sm_tumor_expr_data), threads=infercnv.env$GLOBAL_NUM_THREADS), method=hclust_method)
    
    tumor_subcluster_info$hc = hc
    
    heights = hc$height

    grps <- .partition_by_random_smoothed_trees(tumor_name, tumor_expr_data, hclust_method, p_val, window_size,
                                                max_recursion_depth, min_cluster_size_recurse)

            
    tumor_subcluster_info$subclusters = list()
    
    ordered_idx = tumor_group_idx[hc$order]
    s = split(grps,grps)
    
    flog.info(sprintf("cut tree into: %g groups", length(s)))
    
    start_idx = 1
    for (split_subcluster in names(s)) {
        flog.info(sprintf("-processing %s,%s", tumor_name, split_subcluster))

        split_subcluster_cell_names = names(s[[split_subcluster]])

        if (! all(split_subcluster_cell_names %in% names(tumor_group_idx)) ) {
            stop("Error: .single_tumor_subclustering_smoothed_tree(), not all subcluster cell names were in the tumor group names")
        }
        
        subcluster_indices = tumor_group_idx[ which(names(tumor_group_idx) %in% split_subcluster_cell_names) ]
        
        tumor_subcluster_info$subclusters[[ split_subcluster ]] = subcluster_indices
        
    }
    
    return(tumor_subcluster_info)
}


## Random Trees

.partition_by_random_smoothed_trees <- function(tumor_name, tumor_expr_data, hclust_method, p_val, window_size,
                                                max_recursion_depth, min_cluster_size_recurse) {

    grps <- rep(sprintf("%s.%d", tumor_name, 1), ncol(tumor_expr_data))
    names(grps) <- colnames(tumor_expr_data)

    grps <- .single_tumor_subclustering_recursive_random_smoothed_trees(tumor_expr_data, hclust_method, p_val, grps, window_size,
                                                                        max_recursion_depth, min_cluster_size_recurse)

    
    return(grps)

}


.single_tumor_subclustering_recursive_random_smoothed_trees <- function(tumor_expr_data, hclust_method, p_val, grps.adj, window_size,
                                                                        max_recursion_depth,
                                                                        min_cluster_size_recurse,
                                                                        recursion_depth=1) {


    if (recursion_depth > max_recursion_depth) {
        flog.warn("-not exceeding max recursion depth.")
        return(grps.adj)
    }
    
    tumor_clade_name = unique(grps.adj[names(grps.adj) %in% colnames(tumor_expr_data)])
    message("unique tumor clade name: ", tumor_clade_name)
    if (length(tumor_clade_name) > 1) {
        stop("Error, found too many names in current clade")
    }

    rand_params_info = .parameterize_random_cluster_heights_smoothed_trees(tumor_expr_data, hclust_method, window_size)
    
    h_obs = rand_params_info$h_obs
    h = h_obs$height
    max_height = rand_params_info$max_h
    
    max_height_pval = 1
    if (max_height > 0) {
        ## important... as some clades can be fully collapsed (all identical entries) with zero heights for all
        e = rand_params_info$ecdf
        max_height_pval = 1- e(max_height)
    }

    #message(sprintf("Lengths(h): %s", paste(h, sep=",", collapse=",")))
    #message(sprintf("max_height_pval: %g", max_height_pval))
    
    if (max_height_pval <= p_val) {
        ## keep on cutting.
        cut_height = mean(c(h[length(h)], h[length(h)-1]))
        flog.info(sprintf("cutting at height: %g",  cut_height))
        grps = cutree(h_obs, h=cut_height)
        print(grps)
        uniqgrps = unique(grps)
        
        message("unique grps: ", paste0(uniqgrps, sep=",", collapse=","))

        if (all(sapply(uniqgrps, function(grp) {
            (sum(grps==grp) < min_cluster_size_recurse)
        } ))) {
            flog.warn("none of the split subclusters exceed min cluster size. Not recursing here.")
            return(grps.adj)
        }
                
        for (grp in uniqgrps) {
            grp_idx = which(grps==grp)
            
            message(sprintf("grp: %s  contains idx: %s", grp, paste(grp_idx,sep=",", collapse=","))) 
            df = tumor_expr_data[,grp_idx,drop=FALSE]
            ## define subset.
            subset_cell_names = colnames(df)
            
            subset_clade_name = sprintf("%s.%d", tumor_clade_name, grp)
            message(sprintf("subset_clade_name: %s", subset_clade_name));
            grps.adj[names(grps.adj) %in% subset_cell_names] <- subset_clade_name
            

            if (length(grp_idx) >= min_cluster_size_recurse) {
                ## recurse
                grps.adj <- .single_tumor_subclustering_recursive_random_smoothed_trees(tumor_expr_data=df,
                                                                                        hclust_method=hclust_method,
                                                                                        p_val=p_val,
                                                                                        grps.adj=grps.adj,
                                                                                        window_size=window_size,
                                                                                        max_recursion_depth=max_recursion_depth,
                                                                                        min_cluster_size_recurse=min_cluster_size_recurse,
                                                                                        recursion_depth = recursion_depth + 1 )
            } else {
                flog.warn(sprintf("%s size of %d is too small to recurse on", subset_clade_name, length(grp_idx)))
            }
            
        }
        
    } else {
        message("No cluster pruning: ", tumor_clade_name)
    }
    
    return(grps.adj)
}


.parameterize_random_cluster_heights_smoothed_trees <- function(expr_matrix, hclust_method, window_size, plot=FALSE) {
    
    ## inspired by: https://www.frontiersin.org/articles/10.3389/fgene.2016.00144/full

    sm_expr_data = apply(expr_matrix, 2, caTools::runmean, k=window_size)
    # sm_expr_data = scale(sm_expr_data, center=TRUE, scale=FALSE)
    sm_expr_data = .center_columns(sm_expr_data, 'median')

    
    d = parallelDist(t(sm_expr_data), threads=infercnv.env$GLOBAL_NUM_THREADS)
    
    h_obs = hclust(d, method=hclust_method)
    
        
    # permute by chromosomes
    permute_col_vals <- function(df) {
        ## cells as rows, features as columns
        
        num_cells = nrow(df)

        for (i in seq(ncol(df) ) ) {
            
            df[, i] = df[sample(x=seq_len(num_cells), size=num_cells, replace=FALSE), i]
        }
        
        df
    }
    
    
    flog.info(sprintf("random trees, using %g parallel threads", infercnv.env$GLOBAL_NUM_THREADS))

    available_cores = future::availableCores(which="max")

    if (infercnv.env$GLOBAL_NUM_THREADS > available_cores) {
        flog.warn(sprintf("not enough cores available, setting to num avail cores: %g", available_cores))
        infercnv.env$GLOBAL_NUM_THREADS <- available_cores
    }
    
    # library(doParallel)
    registerDoParallel(cores=infercnv.env$GLOBAL_NUM_THREADS)
    num_rand_iters=100
    max_rand_heights <- foreach (i=seq_len(num_rand_iters)) %dopar% {
        #message("rand iteration: ", i)
        
        rand.tumor.expr.data = t(permute_col_vals( t(expr_matrix) ))
        
        ## smooth it and re-center:
        sm.rand.tumor.expr.data = apply(rand.tumor.expr.data, 2, caTools::runmean, k=window_size)
        # sm.rand.tumor.expr.data = scale(sm.rand.tumor.expr.data, center=TRUE, scale=FALSE)
        sm.rand.tumor.expr.data = .center_columns(sm.rand.tumor.expr.data, 'median')
        
        rand.dist = parallelDist(t(sm.rand.tumor.expr.data), threads=infercnv.env$GLOBAL_NUM_THREADS)
        h_rand <- hclust(rand.dist, method=hclust_method)
        max_rand_height <- max(h_rand$height)

        max_rand_height
    }
    
    max_rand_heights <- as.numeric(max_rand_heights)
    
    h = h_obs$height
    
    max_height = max(h)
    
    message(sprintf("Lengths for original tree branches (h): %s", paste(h, sep=",", collapse=",")))
    message(sprintf("Max height: %g", max_height))

    message(sprintf("Lengths for max heights: %s", paste(max_rand_heights, sep=",", collapse=",")))
    
    e = ecdf(max_rand_heights)
    
    pval = 1- e(max_height)
    message(sprintf("pval: %g", pval))
    
    params_list <- list(h_obs=h_obs,
                        max_h=max_height,
                        rand_max_height_dist=max_rand_heights,
                        ecdf=e
                        )
    
    if (plot) {
        .plot_tree_height_dist(params_list)
    }
    
    
    return(params_list)
    
}


.plot_tree_height_dist <- function(params_list, plot_title='tree_heights') {

    mf = par(mfrow=(c(2,1)))
    
    ## density plot
    rand_height_density = density(params_list$rand_max_height_dist)
    
    xlim=range(params_list$max_h, rand_height_density$x)
    ylim=range(rand_height_density$y)
    plot(rand_height_density, xlim=xlim, ylim=ylim, main=paste(plot_title, "density"))
    abline(v=params_list$max_h, col='red')

        
    ## plot the clustering
    h_obs = params_list$h_obs
    h_obs$labels <- NULL #because they're too long to display
    plot(h_obs)
                
    par(mf)
        
}

.get_tree_height_via_ecdf <- function(p_val, params_list) {
    
    h = quantile(params_list$ecdf, probs=1-p_val)

    return(h)
}



find_DE_stat_significance <- function(normal_matrix, tumor_matrix) {
    
    run_t_test<- function(idx) {
        vals1 = unlist(normal_matrix[idx,,drop=TRUE])
        vals2 = unlist(tumor_matrix[idx,,drop=TRUE])
        
        ## useful way of handling tests that may fail:
        ## https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html

        res = try(t.test(vals1, vals2), silent=TRUE)
        
        if (is(res, "try-error")) return(NA) else return(res$p.value)
        
    }

    pvals = sapply(seq(nrow(normal_matrix)), run_t_test)

    return(pvals)
}


