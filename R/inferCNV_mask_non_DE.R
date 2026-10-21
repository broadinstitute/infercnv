#' @include inferCNV.R
NULL

#' @title mask_non_DE_genes_basic()
#'
#' @description Mask gene expression in infercnv_obj based on expression values that are found
#' statistically significantly different between tumor/normal based on specified
#' statistical test.
#'
#' Uses the expression data stored in infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv object
#'
#' @param p_val_thresh  p-value threshold for assigning DE statistically significant 
#'
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.
#'
#' @param center_val value to assign to those genes that are not found to be statistically DE.
#'
#' @param require_DE_all_normals mask gene if found significantly DE in each normal comparison (default="any") options("any", "most", "all")
#'
#' @return infercnv_obj
#' 
#' @keywords internal
#' @noRd
#'

mask_non_DE_genes_basic <- function(infercnv_obj,
                                    p_val_thresh = 0.05,
                                    test.use="wilcoxon",
                                    center_val=mean(infercnv_obj@expr.data),
                                    require_DE_all_normals="any"
                                    ) {
    
    tumor_groupings = infercnv_obj@observation_grouped_cell_indices

        
    all_DE_results = get_DE_genes_basic(infercnv_obj,
                                        p_val_thresh=p_val_thresh,
                                        test.use=test.use)

    
    #save('all_DE_results', file='all_DE_results.obj')
    
    
    infercnv_obj <- .mask_DE_genes(infercnv_obj,
                                   all_DE_results=all_DE_results,
                                   mask_val=center_val,
                                   require_DE_all_normals=require_DE_all_normals)
    
    return(infercnv_obj)
}




#' @title .mask_DE_genes()
#'
#' @description private function that does the actual masking of expression values in the matrix
#' according to the specified mask value
#'
#' @param infercnv_obj infercnv_object
#'
#' @param all_DE_results  DE results list with structure:
#'                           all_DE_results[[ tumor_type ]] = list(tumor=tumor_type, de_genes=c(geneA, geneB, ...))
#'
#' @param mask_val float value to assign to genes that are in the mask set (complement of the DE gene set)
#'
#' @param min_cluster_size_mask clusters smaller than this size are automatically retained (unmasked). default=5
#' 
#' @return infercnv_obj
#' 
#' @keywords internal
#' @noRd
#'

.mask_DE_genes <- function(infercnv_obj,
                           all_DE_results,
                           mask_val,
                           require_DE_all_normals, # any, most, all
                           min_cluster_size_mask=5) {
        
    all_DE_genes_matrix = matrix(data=0, nrow=nrow(infercnv_obj@expr.data),
                                 ncol=ncol(infercnv_obj@expr.data),
                                 dimnames = list(rownames(infercnv_obj@expr.data),
                                                 colnames(infercnv_obj@expr.data) ) )


    num_normal_types = length(names(infercnv_obj@reference_grouped_cell_indices))

    ## turn on all normal genes
    all_DE_genes_matrix[, unlist(infercnv_obj@reference_grouped_cell_indices)] = num_normal_types

    ## retain small clusters that are unlikely to show up as DE
    #for (tumor_type in names(tumor_groupings)) {
    #    tumor_idx = tumor_groupings[[tumor_type]]
    #    if (length(tumor_idx) < min_cluster_size_mask) {
    #        all_DE_genes_matrix[, tumor_idx] = num_normal_types
    #    }
    #}
    for (DE_results in all_DE_results) {
        if (length(DE_results$tumor_indices) < min_cluster_size_mask) {
            all_DE_genes_matrix[, DE_results$tumor_indices] = num_normal_types
        }
    }
        
    for (DE_results in all_DE_results) {
        # tumor_type = DE_results$tumor
        genes = DE_results$de_genes
        
        gene_idx = rownames(all_DE_genes_matrix) %in% genes
        # cell_idx = tumor_groupings[[ tumor_type ]]
        cell_idx = DE_results$tumor_indices

        if (length(cell_idx) >= min_cluster_size_mask) {
            all_DE_genes_matrix[gene_idx, cell_idx] = all_DE_genes_matrix[gene_idx, cell_idx] + 1
        }
    }
    
    if (require_DE_all_normals == "all") {
        ## must be found in each of the tumor vs (normal_1, normal_2, ..., normal_N) DE comparisons to not be masked.
        infercnv_obj@expr.data[ all_DE_genes_matrix != num_normal_types ] = mask_val
    } else if ( require_DE_all_normals == "most") {
        infercnv_obj@expr.data[ all_DE_genes_matrix < num_normal_types/2 ] = mask_val
    } else if ( require_DE_all_normals == "any") {
        ## masking if not found DE in any comparison
        infercnv_obj@expr.data[ all_DE_genes_matrix == 0 ] = mask_val
    } else {
        stop(sprintf("Error, not recognizing require_DE_all_normals=%s", require_DE_all_normals))
    }
    
    return(infercnv_obj)
    
}



#' @title get_DE_genes_basic
#'
#' @description Retrieves genes identified as significantly DE based on the corresponding statistical test as applied
#' to the infercnv_obj@expr.data with tumor/normal comparison
#'
#' @param infercnv_obj infercnv object
#'
#' @param p_val_thresh maximum p-value for defining statistical signficance
#'
#' @param test.use  statistical test to use. (default: "t" for t-test)
#' 
#' @return all_DE_results[[condition_pair]] = list(tumor=tumor_type,
#'                                                 normal=normal_type,
#'                                                 pvals=pvals,
#'                                                 de_genes=genes)
#'
#' @keywords internal
#' @noRd
#'

get_DE_genes_basic <- function(infercnv_obj,
                               p_val_thresh = 0.05,
                               test.use="wilcoxon" # other options:  (wilcoxon, t, perm)
                               ) {
        
    all_DE_results = list()

    statfxns = list()
    statfxns[[ "t" ]] <- function(x, idx1, idx2) {
        vals1 = x[idx1]
        vals2 = x[idx2]

        # useful way of handling tests that may fail:
        # https://stat.ethz.ch/pipermail/r-help/2008-February/154167.html
        res = try(t.test(vals1, vals2), silent=TRUE)
        
        if (is(res, "try-error")) return(NA) else return(res$p.value)
        
    }

    statfxns[[ "perm" ]] <- function(x, idx1, idx2) {

        vals1 = x[idx1]
        vals2 = x[idx2]

        allvals = c(vals1, vals2)
        facts = factor(rep(c("A","B"), c(length(vals1), length(vals2))))

        perm = coin::oneway_test(allvals ~ facts)
        pval = coin::pvalue(perm)
        return(pval)
    }

    statfxns[[ "wilcoxon" ]] <- function(x, idx1, idx2) {

        vals1 = x[idx1]
        vals2 = x[idx2]

        ## force break ties by adding random noise
        vals1 = vals1 + rnorm(n=length(vals1), mean=0.0001, sd=0.0001)
        vals2 = vals2 + rnorm(n=length(vals2), mean=0.0001, sd=0.0001)
        
        w = wilcox.test(vals1, vals2)

        return(w$p.value)
    }
    
    
    statfxn = statfxns[[ test.use ]]

    ## Find DE genes by comparing the mutant types to normal types
    normal_types = names(infercnv_obj@reference_grouped_cell_indices)
        
    ## turn on only DE genes in tumors
    tumor_groupings = infercnv_obj@observation_grouped_cell_indices
    # observation subclusters are stored per annotation group, or under "all_observations" when
    # cluster_by_groups=FALSE; without subclusters, each annotation group is tested as a whole.
    # (Previously only per-group names were looked up, so with cluster_by_groups=FALSE no tumour cells
    # were tested and all of them were masked; unnamed index vectors were also skipped.)
    subclusters = infercnv_obj@tumor_subclusters[["subclusters"]]
    obs_entries = intersect(names(subclusters), c(names(tumor_groupings), "all_observations"))
    if (length(obs_entries) == 0) {
        subclusters = tumor_groupings
        obs_entries = names(tumor_groupings)
    }
    for (tumor_type in obs_entries) {

        indices = subclusters[[ tumor_type ]]
        if(is.list(indices)) {
            tumor_indices_list = indices
        }
        else { # is.vector(indices)
            tumor_indices_list = list(indices)
            names(tumor_indices_list) = tumor_type
        }
        
        for (tumor_indices_name in names(tumor_indices_list)) {
            tumor_indices = tumor_indices_list[[ tumor_indices_name ]]
#    for (tumor_type in names(tumor_groupings)) {
#
#        tumor_indices = tumor_groupings[[ tumor_type ]] 
        
            for (normal_type in normal_types) {
                flog.info(sprintf("Finding DE genes between %s and %s", tumor_indices_name, normal_type))
    
                normal_indices = infercnv_obj@reference_grouped_cell_indices[[ normal_type ]]
    
    
                pvals = apply(infercnv_obj@expr.data, 1, statfxn, idx1=normal_indices, idx2=tumor_indices)
                pvals = unlist(pvals)
                pvals = p.adjust(pvals, method="BH")
                
                names(pvals) = rownames(infercnv_obj@expr.data)
    
                genes = names(pvals)[pvals<p_val_thresh]
                
                flog.info(sprintf("Found %d genes / %d total as DE", length(genes), length(pvals)))
                
                condition_pair = paste(tumor_indices_name, normal_type, sep=",")
                
                all_DE_results[[condition_pair]] = list(tumor_indices=tumor_indices,
                                                        normal=normal_type,
                                                        pvals=pvals,
                                                        de_genes=genes)
                
                
            }
        }
    }
    
    return(all_DE_results)

}

