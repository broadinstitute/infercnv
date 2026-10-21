
#' @title apply_median_filtering
#' 
#' @description Apply a median filtering to the expression matrix within each tumor bounds
#'
#' @param infercnv_obj infercnv_object
#' 
#' @param window_size Size of the window side centered on the data point to filter (default = 7).
#' 
#' @param on_observations  boolean (default=TRUE), run on observations data (tumor cells).
#'
#' @param on_references  boolean (default=TRUE), run on references (normal cells).
#'
#' @return infercnv_obj with median filtering applied to observations
#'
#' @export
#'
#' @examples
#' # data(infercnv_data_example)
#' # data(infercnv_annots_example)
#' # data(infercnv_genes_example)
#'
#' # infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_example, 
#' #                                                           gene_order_file=infercnv_genes_example,
#' #                                                           annotations_file=infercnv_annots_example,
#' #                                                           ref_group_names=c("normal"))
#'
#' # infercnv_object_example <- infercnv::run(infercnv_object_example,
#' #                                          cutoff=1,
#' #                                          out_dir=tempfile(), 
#' #                                          cluster_by_groups=TRUE, 
#' #                                          denoise=TRUE,
#' #                                          HMM=FALSE,
#' #                                          num_threads=2,
#' #                                          no_plot=TRUE)
#'
#' data(infercnv_object_example)
#'
#' infercnv_object_example <- infercnv::apply_median_filtering(infercnv_object_example)
#' # plot result object
#'

apply_median_filtering <- function(infercnv_obj,
                                   window_size=7,
                                   on_observations=TRUE,
                                   on_references=TRUE) {

    if (window_size%%2 != 1 | window_size < 2) {
      flog.error("::apply_median_filtering: Error, window_size is an even or < 2. Please specify an odd number >= 3.")
    }
    
    half_window = (window_size - 1) / 2
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    if (on_observations) {
        for (tumor_type in names(infercnv_obj@observation_grouped_cell_indices)) {
            
            tumor_indices_list = infercnv_obj@tumor_subclusters[["subclusters"]][[ tumor_type ]]
            # filter neighbouring cells in the order of the subclustering dendrogram (as in the heatmap),
            # not in matrix column order
            hc = infercnv_obj@tumor_subclusters[["hc"]][[ tumor_type ]]
            hc_cell_order = NULL
            if (! is.null(hc) && ! is.null(hc$labels)) {
                hc_cell_order = match(hc$labels[hc$order], colnames(infercnv_obj@expr.data))
            }
            
            for (tumor_indices in tumor_indices_list) {
                if (! is.null(hc_cell_order) && all(tumor_indices %in% hc_cell_order)) {
                    tumor_indices = hc_cell_order[hc_cell_order %in% tumor_indices]
                }
                for (chr in chrs) {
                    chr_genes_indices = which(gene_chr_listing == chr)
                    working_data = infercnv_obj@expr.data[chr_genes_indices, tumor_indices, drop=FALSE]
                    
                    infercnv_obj@expr.data[chr_genes_indices, tumor_indices] = .median_filter(data=working_data,
                                                                                              window_size=window_size,
                                                                                              half_window=half_window)
                }
            }
        }
    }
    
    if (on_references) {
        for (ref_indices in infercnv_obj@reference_grouped_cell_indices) {
            for (chr in chrs) {
                chr_genes_indices = which(gene_chr_listing == chr)
                working_data = infercnv_obj@expr.data[chr_genes_indices, ref_indices, drop=FALSE]
                
                infercnv_obj@expr.data[chr_genes_indices, ref_indices] = .median_filter(data=working_data,
                                                                                        window_size=window_size,
                                                                                        half_window=half_window)
            }
        }
    }
    
    return(infercnv_obj)
}


.median_filter <- function(data,
                           window_size,
                           half_window) {

    xdim = dim(data)[1]
    ydim = dim(data)[2]
    results = data
    
    # if (xdim >= window_size & ydim >= window_size) {
        for (posx in seq_len(xdim)) {
            # window of window_size centred on posx, truncated at the edges (previously half_window + 1 on each side)
            posxa <- max(1, posx - half_window)
            posxb <- min(xdim, posx + half_window)
            for (posy in seq_len(ydim)) {
                posya <- max(1, posy - half_window)
                posyb <- min(ydim, posy + half_window)
                results[posx, posy] = median(data[posxa:posxb, posya:posyb])
            }
        }
    #}

    return(results)
}

