
#' @description Add meta.data about CNAs to a Seurat object from an infercnv_obj
#'
#' @title add_to_seurat()
#'
#' @param seurat_obj Seurat object to add meta.data to (default: NULL)
#'
#' @param infercnv_output_path Path to the output folder of the infercnv run to use
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @param bp_tolerance How many bp of tolerance to have around feature start/end positions for top_n largest CNVs.
#'
#' @return seurat_obj
#'
#' @export
#'

add_to_seurat <- function(seurat_obj = NULL,
                          infercnv_output_path,
                          top_n = 10,
                          bp_tolerance = 2000000) {
    lfiles <- list.files(infercnv_output_path, full.names = FALSE)
    
    if (!file.exists(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))) {
        flog.warn(sprintf("::Could not find \"run.final.infercnv_obj\" file at: %s", paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep)))
        stop()
    }
    infercnv_obj = readRDS(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))
    
    if (is.null(seurat_obj)) {
        flog.info("No Seurat object provided, will only write metadata matrix.")
    }
    else if(!(setequal(row.names(seurat_obj@meta.data), colnames(infercnv_obj@expr.data)) ||
         setequal(colnames(seurat_obj@assays$RNA), colnames(infercnv_obj@expr.data)))) {
        flog.warn("::Cell names in Seurat object and infercnv results do not match")
        stop()
    }

    # all(colnames(infercnv_obj@expr.data)[match(colnames(seurat_obj@assays$RNA), colnames(infercnv_obj@expr.data))] == colnames(seurat_obj@assays$RNA))
    if (is.null(seurat_obj)) {
        cell_ordering_match = seq_len(ncol(infercnv_obj@expr.data))
    }
    else {
        cell_ordering_match = match(colnames(seurat_obj@assays$RNA), colnames(infercnv_obj@expr.data))
    }
    
    ## add check that data row/col names match seurat obj
    analysis_mode_pattern = "samples"
    if (!is.null(infercnv_obj@options$analysis_mode)) {
            analysis_mode_pattern = gsub('[\"]', '', infercnv_obj@options$analysis_mode)   ## gsub to remove _"\"_ around the name of the options if it is given as string and not through a variable
    }

    by_cells = FALSE
    if (!is.null(infercnv_obj@options$HMM_report_by) && infercnv_obj@options$HMM_report_by == "cell") {
        by_cells = TRUE
    }

    scaling_factor = 2
    if (any(grep(lfiles, pattern=paste0("HMM_CNV_predictions.HMM.*", analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
        ###### states used to be 0/0.5/1/1.5/2, they are now 1/2/3/4/5/6
        # scaling_factor = 1
        if (any(grep(lfiles, pattern=paste0("HMM_CNV_predictions.HMMi6.*", analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
            center_state = 3
        }
        else if (any(grep(lfiles, pattern=paste0("HMM_CNV_predictions.HMMi3.*", analysis_mode_pattern, ".Pnorm_0.[0-9]+")))) {
            center_state = 2
        }
        else {
            flog.warn("::Found filtered HMM predictions output, but they do not match any known model type.")
            stop()
        }
        # sort to take lowest BayesProb if there are multiple
        regions = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles, pattern=paste0("HMM_CNV_predictions.HMMi[36].*", analysis_mode_pattern, ".Pnorm_0.[0-9]+.pred_cnv_regions.dat"))])[1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, sort(lfiles[grep(lfiles, pattern=paste0("HMM_CNV_predictions.HMMi[36].*", analysis_mode_pattern, ".Pnorm_0.[0-9]+.pred_cnv_genes.dat"))])[1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        # from_pbayes()
    }
    else if (any(grep(lfiles, pattern = paste0("17_HMM_predHMM.*", analysis_mode_pattern)))) {
        ###### states are 1/2/3/4/5/6
        # scaling_factor = 2
        if (any(grep(lfiles, pattern = paste0("17_HMM_predHMMi6.*", analysis_mode_pattern)))) {
            center_state = 3
        }
        else if (any(grep(lfiles, pattern = paste0("17_HMM_predHMMi3.*", analysis_mode_pattern)))) {
            center_state = 2
        }
        else {
            flog.warn("::Found HMM predictions output, but they do not match any known model type")
            stop()
        }
        regions = read.table(paste(infercnv_output_path, lfiles[grep(lfiles, pattern=paste0("17_HMM_predHMMi[36].*", analysis_mode_pattern, ".pred_cnv_regions.dat"))][1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        hmm_genes = read.table(paste(infercnv_output_path, lfiles[grep(lfiles, pattern=paste0("17_HMM_predHMMi[36].*", analysis_mode_pattern, ".pred_cnv_genes.dat"))][1], sep=.Platform$file.sep), sep="\t", header=TRUE, check.names=FALSE)
        # from_hmm()
    }
    else {
        flog.warn(sprintf("::Could not find any HMM predictions outputs at: %s", infercnv_output_path))
        flog.warn("::add_to_seurat works by transcribing (and copying if a Seurat object is provided) the HMM results from infercnv, so it cannot be used without having run the HMM.")
        stop()
    }
    
    features_to_add <- .get_features(infercnv_obj = infercnv_obj,
                                     infercnv_output_path = infercnv_output_path,
                                     regions = regions, 
                                     hmm_genes = hmm_genes, 
                                     center_state = center_state, 
                                     scaling_factor = scaling_factor,
                                     by_cells = by_cells,
                                     top_n = top_n,
                                     bp_tolerance = bp_tolerance)
    if (!is.null(seurat_obj)) {
        for (lv in levels(infercnv_obj@gene_order$chr)) {
            seurat_obj@meta.data[[paste0("has_cnv_", lv)]] = features_to_add$feature_vector_chrs_has_cnv[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("has_loss_", lv)]] = features_to_add$feature_vector_chrs_has_loss[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("has_dupli_", lv)]] = features_to_add$feature_vector_chrs_has_dupli[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_cnv_", lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_loss_", lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_dupli_", lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_scaled_cnv_", lv)]] = features_to_add$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_scaled_loss_", lv)]] = features_to_add$feature_vector_chrs_gene_loss_proportion_scaled[[lv]][cell_ordering_match][cell_ordering_match]
            seurat_obj@meta.data[[paste0("proportion_scaled_dupli_", lv)]] = features_to_add$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]][cell_ordering_match]
        }
        
        for (n in names(features_to_add)[grep(names(features_to_add), pattern = "top_")] ) {
            seurat_obj@meta.data[[n]] = features_to_add[[n]][cell_ordering_match]
        }
    }
    
    out_mat = matrix(NA, ncol=((9 * length(levels(infercnv_obj@gene_order$chr))) + length(features_to_add) - 9), nrow=ncol(infercnv_obj@expr.data))
    out_mat_feature_names = vector("character", ((9 * length(levels(infercnv_obj@gene_order$chr))) + length(features_to_add) - 9))
    
    i = 1
    for (lv in levels(infercnv_obj@gene_order$chr)) {
        
        out_mat[, i] = features_to_add$feature_vector_chrs_has_cnv[[lv]]
        out_mat[, i+1] = features_to_add$feature_vector_chrs_has_loss[[lv]]
        out_mat[, i+2] = features_to_add$feature_vector_chrs_has_dupli[[lv]]
        out_mat[, i+3] = features_to_add$feature_vector_chrs_gene_cnv_proportion[[lv]]
        out_mat[, i+4] = features_to_add$feature_vector_chrs_gene_loss_proportion[[lv]]
        out_mat[, i+5] = features_to_add$feature_vector_chrs_gene_dupli_proportion[[lv]]
        out_mat[, i+6] = features_to_add$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]]
        out_mat[, i+7] = features_to_add$feature_vector_chrs_gene_loss_proportion_scaled[[lv]]
        out_mat[, i+8] = features_to_add$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]]
        
        out_mat_feature_names[i:(i+8)] = c(paste0("has_cnv_", lv), paste0("has_loss_", lv), paste0("has_dupli_", lv), paste0("proportion_cnv_", lv), paste0("proportion_loss_", lv), paste0("proportion_dupli_", lv), paste0("proportion_scaled_cnv_", lv), paste0("proportion_scaled_loss_", lv), paste0("proportion_scaled_dupli_", lv))
        
        i = i + 9
    }
    
    for (n in names(features_to_add)[grep(names(features_to_add), pattern = "top_")] ) {
        out_mat[, i] = features_to_add[[n]]
        out_mat_feature_names[i] = n
        i = i + 1
    }
    
    colnames(out_mat) = out_mat_feature_names
    row.names(out_mat) = colnames(infercnv_obj@expr.data)
    
    write.table(out_mat, paste(infercnv_output_path, "map_metadata_from_infercnv.txt", sep=.Platform$file.sep) , quote=FALSE, sep="\t")
    
    return(seurat_obj)
}


#' @title .get_features()
#'
#' @description Get data from infercnv objects to add to Seurat meta.data
#'
#' @param infercnv_obj infercnv hmm object
#'
#' @param infercnv_output_path infercnv output folder
#'
#' @param regions Table with predicted CNAs regions from the HMM model
#'
#' @param hmm_genes Table with predicted CNAs genes from the HMM model
#'
#' @param center_state Value that represents the neutral state in the HMM results.
#'
#' @param scaling_factor Factor to multiply divergence from center_state to get CNA amplitude.
#'
#' @param by_cells Whether the HMM predictions are given by cells or not.
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @return all_features A list of all the calculated meta.data to add.
#'
#' @keywords internal
#' @noRd
#'
.get_features <- function(infercnv_obj, infercnv_output_path, regions, hmm_genes, center_state, scaling_factor, by_cells, top_n, bp_tolerance) {
    
    chr_gene_count = table(infercnv_obj@gene_order$chr)
    
    # features templates for initialization
    double_feature_vector = vector(mode="double", length=ncol(infercnv_obj@expr.data))
    names(double_feature_vector) = colnames(infercnv_obj@expr.data)
    logical_feature_vector = vector(mode="logical", length=ncol(infercnv_obj@expr.data))
    names(logical_feature_vector) = colnames(infercnv_obj@expr.data)
    
    # initialize features lists
    all_features = c()
    all_features$feature_vector_chrs_has_cnv = c()
    all_features$feature_vector_chrs_has_dupli = c()
    all_features$feature_vector_chrs_has_loss = c()
    all_features$feature_vector_chrs_gene_cnv_proportion = c()
    all_features$feature_vector_chrs_gene_dupli_proportion = c()
    all_features$feature_vector_chrs_gene_loss_proportion = c()
    all_features$feature_vector_chrs_gene_cnv_proportion_scaled = c()
    all_features$feature_vector_chrs_gene_dupli_proportion_scaled = c()
    all_features$feature_vector_chrs_gene_loss_proportion_scaled = c()
    for (lv in levels(infercnv_obj@gene_order$chr)) {
        all_features$feature_vector_chrs_has_cnv[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_has_dupli[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_has_loss[[lv]] = logical_feature_vector
        all_features$feature_vector_chrs_gene_cnv_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_dupli_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_loss_proportion[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_cnv_proportion_scaled[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_dupli_proportion_scaled[[lv]] = double_feature_vector
        all_features$feature_vector_chrs_gene_loss_proportion_scaled[[lv]] = double_feature_vector
    }
    
    # map for top_n mapping
    subclust_name_to_clust = list()

    if (!by_cells) {
        for (clust in names(infercnv_obj@tumor_subclusters$subclusters)) {
            for (subclust in names(infercnv_obj@tumor_subclusters$subclusters[[clust]])) {
                subclust_name = paste(clust, subclust, sep=".")
                subclust_name_to_clust[[subclust_name]] = c(clust, subclust)
                res = regions[regions$cell_group_name == subclust_name, , drop=FALSE]
                gres = hmm_genes[hmm_genes$cell_group_name == subclust_name, , drop=FALSE]
                if (nrow(res) > 0) {
                    for (c in unique(res$chr)) {
                        all_features$feature_vector_chrs_has_cnv[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                        all_features$feature_vector_chrs_gene_cnv_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(gres$chr == c)) / chr_gene_count[[c]])
                        all_features$feature_vector_chrs_gene_cnv_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (sum(abs(gres[(which(gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                    }
                    
                    sub_gres = gres[gres$state < center_state, ]
                    for (c in unique(sub_gres$chr)) {
                        all_features$feature_vector_chrs_has_loss[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                        all_features$feature_vector_chrs_gene_loss_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                        all_features$feature_vector_chrs_gene_loss_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (abs(sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                    }
                    
                    sub_gres = gres[gres$state > center_state, ]
                    for (c in unique(sub_gres$chr)) {
                        all_features$feature_vector_chrs_has_dupli[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                        all_features$feature_vector_chrs_gene_dupli_proportion[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                        all_features$feature_vector_chrs_gene_dupli_proportion_scaled[[c]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = (sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state) / (chr_gene_count[[c]] * scaling_factor))
                    }
                }
            }
        }
    } else {
        for (cell_id in seq_len(ncol(infercnv_obj@expr.data))) {
            cell_name = colnames(infercnv_obj@expr.data)[cell_id]
            res = regions[regions$cell_group_name == cell_name, , drop=FALSE]
            gres = hmm_genes[hmm_genes$cell_group_name == cell_name, , drop=FALSE]
            if (nrow(res) > 0) {
                for (c in unique(res$chr)) {
                    all_features$feature_vector_chrs_has_cnv[[c]][cell_name] = TRUE
                    all_features$feature_vector_chrs_gene_cnv_proportion[[c]][cell_name] = (length(which(gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_cnv_proportion_scaled[[c]][cell_name] = (sum(abs(gres[(which(gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                }
                
                sub_gres = gres[gres$state < center_state, ]
                for (c in unique(sub_gres$chr)) {
                    all_features$feature_vector_chrs_has_loss[[c]][cell_name] = TRUE
                    all_features$feature_vector_chrs_gene_loss_proportion[[c]][cell_name] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_loss_proportion_scaled[[c]][cell_name] = (abs(sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state)) / (chr_gene_count[[c]] * scaling_factor))
                }
                
                sub_gres = gres[gres$state > center_state, ]
                for (c in unique(sub_gres$chr)) {
                    all_features$feature_vector_chrs_has_dupli[[c]][cell_name] = TRUE
                    all_features$feature_vector_chrs_gene_dupli_proportion[[c]][cell_name] = (length(which(sub_gres$chr == c)) / chr_gene_count[[c]])
                    all_features$feature_vector_chrs_gene_dupli_proportion_scaled[[c]][cell_name] = (sum(sub_gres[(which(sub_gres$chr == c)), "state"] - center_state) / (chr_gene_count[[c]] * scaling_factor))
                }
            }
        }
    }
    
    # sorted_regions = sort(table(hmm_genes$gene_region_name), decreasing=TRUE)
    sorted_regions_loss = sort(table(hmm_genes$gene_region_name[hmm_genes$state < center_state]), decreasing=TRUE)
    sorted_regions_dupli = sort(table(hmm_genes$gene_region_name[hmm_genes$state > center_state]), decreasing=TRUE)
    
    # top_n_cnv = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions = sorted_regions, top_n = top_n, bp_tolerance = bp_tolerance)
    top_n_loss = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions = sorted_regions_loss, top_n = top_n, bp_tolerance = bp_tolerance)
    top_n_dupli = .get_top_n_regions(hmm_genes = hmm_genes, sorted_regions = sorted_regions_dupli, top_n = top_n, bp_tolerance = bp_tolerance)
    
    # for (i in seq_along(top_n_cnv)) {
    #     feature_name = paste0("top_cnv_", i)
    #     all_features[[feature_name]] = logical_feature_vector
        
    #     for (subclust_name in top_n_cnv[[i]]$subclust_name) {
    #         clust = subclust_name_to_clust[[subclust_name]][1]
    #         subclust = subclust_name_to_clust[[subclust_name]][2]
    #         all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
    #     }
    # }

    to_write = c()
    for (i in seq_along(top_n_loss)) {
        feature_name = paste0("top_loss_", i)
        all_features[[feature_name]] = logical_feature_vector
        line_to_write = c()
        
        if (!by_cells) {
            for (subclust_name in top_n_loss[[i]]$subclust_name) {
                clust = subclust_name_to_clust[[subclust_name]][1]
                subclust = subclust_name_to_clust[[subclust_name]][2]
                all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                # line_to_write = c(line_to_write, names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]]))
                line_to_write = c(line_to_write, paste(subclust_name, names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]]), sep=";"))
            }
        } else {
            for (cell_id in top_n_loss[[i]]$subclust_name) {
                all_features[[feature_name]][cell_id] = TRUE
                line_to_write = c(line_to_write, paste(cell_id, cell_id, sep=";"))
            }
        }

        # to_write = c(to_write, paste(feature_name, paste(top_n_loss[[i]]$subclust_name, sep=" "), paste(line_to_write, sep=" "), sep=";"))
        to_write = c(to_write, paste(feature_name, line_to_write, sep=";"))

    }

    if (is.null(to_write)) {
        to_write = ""
    }

    fileConn<-file(paste(infercnv_output_path, "top_losses.txt", sep=.Platform$file.sep), open="w")
    writeLines(to_write, con=fileConn)
    close(fileConn)

    # to_write = to_write = vector(mode="vector", length=length(top_n_dupli))
    to_write = c()
    for (i in seq_along(top_n_dupli)) {
        feature_name = paste0("top_dupli_", i)
        all_features[[feature_name]] = logical_feature_vector
        line_to_write = c()

        if (!by_cells) {
            for (subclust_name in top_n_dupli[[i]]$subclust_name) {
                clust = subclust_name_to_clust[[subclust_name]][1]
                subclust = subclust_name_to_clust[[subclust_name]][2]
                all_features[[feature_name]][names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]])] = TRUE
                # line_to_write = c(line_to_write, paste(feature_name, paste(top_n_dupli[[i]]$subclust_name, sep=" "), paste(line_to_write, sep=" "), sep=";"))
                 # line_to_write = c(line_to_write, names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]]))
                 line_to_write = c(line_to_write, paste(subclust_name, names(infercnv_obj@tumor_subclusters$subclusters[[clust]][[subclust]]), sep=";"))
            }
        } else {
            for (cell_id in top_n_dupli[[i]]$subclust_name) {
                all_features[[feature_name]][cell_id] = TRUE
                line_to_write = c(line_to_write, paste(cell_id, cell_id, sep=";"))
            }
        }

        # to_write = c(to_write, paste(feature_name, paste(top_n_dupli[i]]$subclust_name, sep=" "), paste(line_to_write, sep=" "), sep=";"))
        to_write = c(to_write, paste(feature_name, line_to_write, sep=";"))

    }

    if (is.null(to_write)) {
        to_write = ""
    }

    fileConn<-file(paste(infercnv_output_path, "top_duplis.txt", sep=.Platform$file.sep), open="w")
    writeLines(to_write, con=fileConn)
    close(fileConn)
    
    return(all_features)
}


#' @title .get_top_n_regions()
#'
#' @description Get top n largest CNA regions in number of genes
#'
#' @param hmm_genes Table with predicted CNAs genes from the HMM model
#'
#' @param sorted_region List of regions sorted by size in number of genes for the CNA type desired (gain/loss/both)
#'
#' @param top_n How many of the largest CNA (in number of genes) to get.
#'
#' @return all_features A list of all the calculated meta.data to add.
#'
#' @keywords internal
#' @noRd
#'
.get_top_n_regions <- function(hmm_genes, sorted_regions, top_n, bp_tolerance) {
    
    if (is.null(nrow(sorted_regions)) || (length(sorted_regions) < 1)) {
        return(c())
    }

    j = 1
    top_regions = vector("list", top_n)
    used_regions = c()
    # flog.debug("sorted regions are:")
    # for(sr in names(sorted_regions)) {
    #     flog.debug(paste(sr, sorted_regions[sr]))
    # }

    for (i in seq_len(nrow(sorted_regions))) {

        if (names(sorted_regions[i]) %in% used_regions) {
            next
        }

        genes_in_region = hmm_genes[which(hmm_genes$gene_region_name %in% names(sorted_regions[i])), ]
        region_chr = genes_in_region$chr[1]
        region_start_low = min(genes_in_region$start)
        region_start_high = region_start_low
        region_end_low = max(genes_in_region$end)
        region_end_high = region_end_low

        to_ignore = which(hmm_genes$gene_region_name %in% used_regions)
        
        if (length(to_ignore) > 0) {
            same_chr = setdiff(which(hmm_genes$chr == region_chr), to_ignore)
        }
        else {
            same_chr = which(hmm_genes$chr == region_chr)
        }

        initial_close = list()
        repeat {

            close_start = same_chr[which((hmm_genes$start[same_chr] <= region_start_high + bp_tolerance) & (hmm_genes$start[same_chr] >= region_start_low - bp_tolerance))]
            close_end = same_chr[which((hmm_genes$end[same_chr] <= region_end_high + bp_tolerance) & (hmm_genes$end[same_chr] >= region_end_low - bp_tolerance))]
            close_start_end = intersect(unique(hmm_genes$gene_region_name[close_start]), unique(hmm_genes$gene_region_name[close_end]))
        
            if ((length(setdiff(close_start_end, initial_close)) == 0) && (length(setdiff(initial_close, close_start_end)) == 0)) {
                break
            }
            else {
                initial_close = close_start_end
                starts = c()
                ends = c()
                for (regi in close_start_end) {
                    starts = c(starts, min(hmm_genes$start[which(hmm_genes$gene_region_name == regi)]))
                    ends = c(ends, max(hmm_genes$end[which(hmm_genes$gene_region_name == regi)]))
                }

                region_start_low = min(starts)
                region_start_high = max(starts)
                region_end_low = min(ends)
                region_end_high = max(ends)
            }
        }

        if (length(close_start_end) > 0) {
            top_regions[[j]]$subclust_names = unique(hmm_genes$cell_group_name[which(hmm_genes$gene_region_name %in% close_start_end)])
            top_regions[[j]]$regions_names = close_start_end
            flog.debug(paste0("top cnv ", j, " is composed of subclusts: "))#, paste(close_start_end, sep="   ")))
            flog.debug(paste(top_regions[[j]]$subclust_names, sep="  "))
            flog.debug("and region names: ")
            flog.debug(paste(top_regions[[j]]$regions_names, sep="  "))
        }
        else {
            flog.error("Did not even find itself, error.")
            stop()
        }

        used_regions = c(used_regions, close_start_end)

        if (length(used_regions) != length(unique(used_regions))) {
            flog.error("Used the same region twice")
            stop()
        }

        j = j + 1

        if (j == (top_n + 1)) {
            break
        }
    }

    return(top_regions[1:(j-1)])
}


# .get_top_n_regions <- function(hmm_genes, sorted_regions, top_n, bp_tolerance) {
#     j = 1
#     previous_region_chr = -1
#     previous_region_start = -1
#     previous_region_end = -1
#     top_regions = vector("list", top_n)
    
#     for (i in seq_len(nrow(sorted_regions))) {
#         genes_in_region = hmm_genes[which(hmm_genes$gene_region_name %in% names(sorted_regions[i])), ]
#         region_chr = genes_in_region$chr[1]
#         region_start = min(genes_in_region$start)
#         region_end = max(genes_in_region$end)
#         # check if the current region is the same as the previous one for a different subcluster or not
#         # if it is, extend the previous assignment without increasing the count of found top hits
#         if (region_chr == previous_region_chr && region_start <= previous_region_start + bp_tolerance && region_start >= previous_region_start - bp_tolerance && region_end <= previous_region_end + bp_tolerance && region_end >= previous_region_end - bp_tolerance) {
#             top_regions[[j]]$subclust_names = c(top_regions[[j]]$subclust_names, genes_in_region$cell_group_name[1])
#             top_regions[[j]]$regions_names = c(top_regions[[j]]$regions_names, genes_in_region$gene_region_name[1])
#         }
#         else {
#             top_regions[[j]]$subclust_names = genes_in_region$cell_group_name[1]
#             top_regions[[j]]$regions_names = genes_in_region$gene_region_name[1]
#             previous_region_chr = region_chr
#             previous_region_start = region_start
#             previous_region_end = region_end
#             j = j + 1
#         }
#         if (j == top_n + 1) {
#             break
#         }
#     }
    
#     if (j < top_n + 1) { # if less non unique regions than top_n
#         top_regions = top_regions[1:j]
#     }
    
#     return(top_regions)
# }




##' @keywords internal
##' @noRd
##'
#make_seurat_from_infercnv_obj <- function(infercnv_obj) {
#    return(CreateSeuratObject(counts = infercnv_obj@count.data, project="infercnv", min.cells = 3, min.features = 200))
#}
#
##' @keywords internal
##' @noRd
##'
#make_seurat_from_infercnv <- function(infercnv_output_path) {
#    if (file.exists(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))) {
#        return(make_seurat_from_infercnv_obj(readRDS(paste(infercnv_output_path, "run.final.infercnv_obj", sep=.Platform$file.sep))))
#    }
#    else {
#        stop()
#    }
#}


