
#' @title sample_object
#' 
#' @description Apply sampling on an infercnv object to reduce the number of cells in it 
#' and allow faster plotting or have all groups take up the same height on the heatmap
#'
#' @param infercnv_obj infercnv_object
#' 
#' @param n_cells Number of cells that should be sampled per group (default = 100).
#' 
#' @param every_n Sample 1 cell every_n cells for each group. If subclusters are defined, 
#' this will make sure that at least one cell per subcluster is sampled. 
#' Requires above_m to be set to work, overriding n_cells parameter.
#' 
#' @param above_m Sample groups that have at least above_m cells. 
#' Requires every_n to be set to work, overriding n_cells parameter
#' 
#' @param on_observations  boolean (default=TRUE), sample observations data (tumor cells).
#'
#' @param on_references  boolean (default=TRUE), sample references (normal cells).
#'
#' @return sampled infercnv_obj 
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
#' infercnv_object_example <- infercnv::sample_object(infercnv_object_example, n_cells=5)
#' # plot result object
#'


sample_object <- function(infercnv_obj,
    n_cells=100,
    every_n=NULL,
    above_m=NULL,
    on_references=TRUE,
    on_observations=TRUE) {

    do_every_n = FALSE
    if (!is.null(every_n) && !is.null(above_m)) {
        if (every_n < 2) {
            flog.error("every_n needs to be at least 2, otherwise nothing will be done.")
            stop("every_n < 2")
        }
        else if (floor(every_n) != every_n) {
            flog.error("every_n needs to be an integer.")
            stop("every_n not an integer")
        }
        else if(!("numeric" %in% is(above_m))) {
            flog.error("above_m needs to be numeric.")
            stop("above_m not numeric")
        }
        do_every_n = TRUE
        n_cells = NULL # mostly to make sure it's not used by inadvertence
    }
    else if (!is.null(every_n) || !is.null(above_m)) {
        flog.info("To use object sampling with every_n and above_m options, please set both. Checking if n_cells is set.")
        # stop("every_n and above_m both need to be set or neither.")
    }
    if (!do_every_n) {
        if (is.null(n_cells) || n_cells < 1) {
            flog.error("Please provide a valid number of cells to sample to.")
            stop("invalid n_cells")
        }
    }

    new_obj <- new(
        Class = "infercnv",
        expr.data = matrix(), 
        count.data = matrix(),
        gene_order = infercnv_obj@gene_order,
        reference_grouped_cell_indices = list(),
        observation_grouped_cell_indices = list(),
        tumor_subclusters = infercnv_obj@tumor_subclusters,  # mainly copied for structure, will get updated
        .hspike = infercnv_obj@.hspike)

    col1 = length(unlist(infercnv_obj@reference_grouped_cell_indices))
    col2 = length(unlist(infercnv_obj@observation_grouped_cell_indices))

    if (do_every_n) {
        if (on_references == TRUE) {
            col1 = sum(vapply(names(infercnv_obj@reference_grouped_cell_indices), function(sample_name) {
                if (length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) > above_m) {
                    as.integer(ceiling(length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) / every_n))
                }
                else {
                    length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])
                }
            }, integer(1)))
        }
        if (on_observations == TRUE) {
            col2 = sum(vapply(names(infercnv_obj@observation_grouped_cell_indices), function(sample_name) {
                if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) > above_m) {
                    as.integer(ceiling(length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) / every_n))
                }
                else {
                    length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
                }
            }, integer(1)))
        }
        # make the result matrix slightly bigger than expected in case of subclusters not getting sampled 
        # this allows to add a cell from each non represented subcluster without having to extend the matrix
        # we can then remove the empty columns from this extension based on i at the end (should never be more than this extension)
        # which will be more efficient than extending the matrix multiple times
        col2 = col2 + length(unlist(infercnv_obj@tumor_subclusters$subclusters, recursive=FALSE)) - length(infercnv_obj@tumor_subclusters$subclusters)
    }
    else {
        if (on_references == TRUE) {
            col1 = length(infercnv_obj@reference_grouped_cell_indices) * n_cells
        }
        if (on_observations == TRUE) {
            col2 = length(infercnv_obj@observation_grouped_cell_indices) * n_cells
        }
    }
    new_obj@expr.data = matrix(nrow=nrow(infercnv_obj@expr.data), ncol=(col1 + col2))
    row.names(new_obj@expr.data) = row.names(infercnv_obj@expr.data)
    futur_colnames = rep("", ncol(new_obj@expr.data))
    new_obj@tumor_subclusters$subclusters=list()

    i = 1
    if (on_references == TRUE) {  ## may need to fix @tumor_subclusters$hc of references at some point

        for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {
            if ((!is.null(n_cells) && length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) >= n_cells)
                 || do_every_n) {  # downsample

                flog.info(paste("Downsampling ", sample_name, sep=""))

                if (!do_every_n) {
                    sampled_indices = sample(infercnv_obj@reference_grouped_cell_indices[[sample_name]], size=n_cells, replace=FALSE)
                }
                else if (length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) > above_m) {  # select 1 every_n because above_m
                    seq_indices = seq(1, length(unlist(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])), every_n)

                    sampled_labels = (infercnv_obj@tumor_subclusters$hc[[sample_name]]$labels[infercnv_obj@tumor_subclusters$hc[[sample_name]]$order])[seq_indices]
                    sampled_indices = match(sampled_labels, colnames(infercnv_obj@expr.data))

                    # check that all subclusters are represented
                    for (subcluster_id in names(infercnv_obj@tumor_subclusters$subcluster[[sample_name]])) {
                        if (!any(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]] %in% sampled_indices)) {
                            sampled_indices = c(sampled_indices, infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]][1])
                        }
                    }
                }
                else { # keep everything because we would select 1 every_n but we are not above_m
                    sampled_indices = infercnv_obj@reference_grouped_cell_indices[[sample_name]]
                }
                futur_colnames[(i:(i + length(sampled_indices) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE])
            }
            else if (!is.null(n_cells)) {  # upsample

                flog.info(paste("Upsampling ", sample_name, sep=""))

                n_copies = floor(n_cells / length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]))
                to_sample = n_cells %% length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])

                pre_sampled_indices = sample(seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]]), size=to_sample, replace=FALSE)
                sampled_indices = sort(c(pre_sampled_indices, rep(seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]]), n_copies)))
                sampled_indices = infercnv_obj@reference_grouped_cell_indices[[sample_name]][sampled_indices]
                ### add check in case there is no hclust for references because it was not set in run() options
                if (!is.null(infercnv_obj@tumor_subclusters$hc[[sample_name]])) {
                    new_data_order = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                            if (x %in% pre_sampled_indices) {
                                rep(x, (n_copies + 1))
                            }
                            else {
                                rep(x, n_copies)
                            }
                        }))

                    str_newick_to_alter = write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[sample_name]]))
                    labels = colnames(infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE])
                    for (k in seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]])) {
                        if (k %in% pre_sampled_indices) {
                            # need n_copies + 1 duplicates, so need n_copies added
                            current_label = labels[k]
                            to_replace = current_label
                            for (l in seq_len(n_copies)) {
                                replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                                str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                                to_replace = paste(current_label, "_", (l+1), sep="")
                            }
                        }
                        else if (n_copies > 1) {
                            # make n_copies duplicates, so need n_copies - 1 added
                            current_label = labels[k]
                            to_replace = current_label
                            for (l in seq_len(n_copies - 1)) {
                                replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                                str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                                to_replace = paste(current_label, "_", (l+1), sep="")
                            }
                        }
                    }
                    new_obj@tumor_subclusters$hc[[sample_name]] = as.hclust(read.tree(text=str_newick_to_alter))
                }
                else {
                    if (length(infercnv_obj@tumor_subclusters$subclusters[[sample_name]]) > 1) {
                        flog.error(paste("No hclust object available for ", sample_name, " even though there is more than 1 cell", sep=""))
                        stop("Missing clustering information")
                    }

                    ### add check that n_copies is at least 2
                    # although if there is not hclust, there should only be 1 cell, so n_copies is at least 2 in that case, otherwise there would be no upsampling
                    new_data_order = rep(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[1]], n_cells)
                    new_suffixes = seq_len(n_cells)

                    current_label = colnames(infercnv_obj@expr.data[, infercnv_obj@tumor_subclusters$subclusters[[sample_name]], drop=FALSE])
                    # str_newick_to_alter = paste(paste(lapply(seq_len((n_copies) - 1), function(l) {
                    #         paste("(", current_label, "_", l, ":0", sep="")
                    #     }), collapse=","), ",", current_label, "_", n_copies, ":0", paste(rep("):0", (n_copies - 2)), collapse=""), ");", sep="")
                    new_obj@tumor_subclusters$hc[[sample_name]] <- list()  # initialize empty object
                    new_obj@tumor_subclusters$hc[[sample_name]]$merge <- matrix(c(-1, -2), ncol=2, byrow=TRUE)
                    for (k in seq_len(n_copies - 2)) {
                        new_obj@tumor_subclusters$hc[[sample_name]]$merge <- rbind(new_obj@tumor_subclusters$hc[[sample_name]]$merge, c(k, -(k+2)))
                    }
                    new_obj@tumor_subclusters$hc[[sample_name]]$height <- rep(1, (n_copies - 1))    # define merge heights
                    new_obj@tumor_subclusters$hc[[sample_name]]$order <- seq_len(n_copies)              # order of leaves(trivial if hand-entered)
                    new_obj@tumor_subclusters$hc[[sample_name]]$labels <- paste(current_label, seq_len(n_copies), sep="_")    # lnew_obj@tumor_subclusters$hc[[sample_name]]bels of lenew_obj@tumor_subclusters$hc[[sample_name]]ves
                    class(new_obj@tumor_subclusters$hc[[sample_name]]) <- "hclust"
                }

                futur_colnames[(i:(i + length(sampled_indices) - 1))] = paste(colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE]), seq_along(sampled_indices), sep="_")  # paste with seq_along to ensure unique labels
            }
            new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = c(i:(i + length(sampled_indices) - 1))
            new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, sampled_indices, drop=FALSE]
            new_obj@reference_grouped_cell_indices[[sample_name]]=c(i:(i + length(sampled_indices) - 1))
            i = i + length(sampled_indices)
        }
    }
    else { # do nothing on references
        for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {
            new_obj@expr.data[, (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE]
            new_obj@reference_grouped_cell_indices[[sample_name]] = (i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))
            new_obj@tumor_subclusters$hc[[sample_name]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]
            new_obj@tumor_subclusters$subclusters[[sample_name]] = infercnv_obj@tumor_subclusters$subclusters[[sample_name]]
            futur_colnames[(i:(i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE])
            i = i + length(infercnv_obj@reference_grouped_cell_indices[[sample_name]])
        }
    }

    if (on_observations == TRUE) {
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {
            if ((!is.null(n_cells) && length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) >= n_cells) || do_every_n) { ## downsample 

                flog.info(paste("Downsampling ", sample_name, sep=""))

                if (!do_every_n) {  # basic random sampling of n_cells
                    sampled_indices = sample(infercnv_obj@observation_grouped_cell_indices[[sample_name]], size=n_cells, replace=FALSE)
                }
                else if (length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) > above_m) {  # select 1 every_n because above_m
                    ## sort(unlist(infercnv_obj_subsample@tumor_subclusters$subclusters[[sample_name]], use.names = FALSE)) == infercnv_obj_subsample@reference_grouped_cell_indices[[sample_name]]
                    seq_indices = seq(1, length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), every_n)
                    sampled_labels = (infercnv_obj@tumor_subclusters$hc[[sample_name]]$labels[infercnv_obj@tumor_subclusters$hc[[sample_name]]$order])[seq_indices]
                    sampled_indices = match(sampled_labels, colnames(infercnv_obj@expr.data))

                    # check that all subclusters are represented
                    for (subcluster_id in names(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])) {
                        if (!any(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]] %in% sampled_indices)) {
                            sampled_indices = c(sampled_indices, infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]][1])
                        }
                    }
                }
                else { # keep everything because we would select 1 every_n but we are not above_m
                    sampled_indices = infercnv_obj@observation_grouped_cell_indices[[sample_name]]
                }
                ## prune tree
                to_prune = colnames(infercnv_obj@expr.data)[setdiff(infercnv_obj@observation_grouped_cell_indices[[sample_name]], sampled_indices)]
                
                flog.info(paste("Pruning ", length(to_prune), " leaves from hclust for ", sample_name, sep=""))

                if (length(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order) > (length(to_prune) + 1) ) {  # more than 2 leaves left, so hclust can exist
                    # new_obj@tumor_subclusters$hc[[sample_name]] = prune(infercnv_obj@tumor_subclusters$hc[[sample_name]], to_prune)
                    new_obj@tumor_subclusters$hc[[sample_name]] = as.hclust(drop.tip(as.phylo(infercnv_obj@tumor_subclusters$hc[[sample_name]]), to_prune))
                }
                else { # 1 or no leaves left, can't have such an hclust object
                    flog.info(paste("Not enough leaves left for hclust to exist for ", sample_name, sep=""))
                    new_obj@tumor_subclusters$hc[[sample_name]] = NULL
                }
                new_obj@expr.data[, (i:(i + length(sampled_indices) - 1))] = infercnv_obj@expr.data[, sampled_indices, drop=FALSE]
                futur_colnames[(i:(i + length(sampled_indices) - 1))] = colnames(infercnv_obj@expr.data[, sampled_indices, drop=FALSE])
            }

            else if (!is.null(n_cells)) { ## upsample

                flog.info(paste("Upsampling ", sample_name, sep=""))

                n_copies = floor(n_cells / length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]))
                to_sample = n_cells %% length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
                pre_sampled_indices = sample(seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]]), size=to_sample, replace=FALSE)

                if (!is.null(infercnv_obj@tumor_subclusters$hc[[sample_name]])) {

                    flog.info(paste("Adding leaves to hclust for ", sample_name, sep=""))

                    new_data_order = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                            if (x %in% pre_sampled_indices) {
                                rep(x, (n_copies + 1))
                            }
                            else {
                                rep(x, n_copies)
                            }
                        }))

                    new_suffixes = unlist(lapply(infercnv_obj@tumor_subclusters$hc[[sample_name]]$order, function(x) {
                            if (x %in% pre_sampled_indices) {
                                seq_len(n_copies + 1)
                            }
                            else {
                                seq_len(n_copies)
                            }
                        }))

                    str_newick_to_alter = write.tree(as.phylo(infercnv_obj@tumor_subclusters$hc[[sample_name]]))

                    # write.tree(base_newick)
                    # plot(read.tree(text=write.tree(base_newick)))
                    # gsub(to_replace, replacement, str_newick_to_alter)

                    labels = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE])
                    for (k in seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]])) {
                        if (k %in% pre_sampled_indices) {
                            # need n_copies + 1 duplicates, so need n_copies added
                            current_label = labels[k]
                            to_replace = current_label
                            for (l in seq_len(n_copies)) {
                                replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                                str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                                to_replace = paste(current_label, "_", (l+1), sep="")
                            }
                        }
                        else if (n_copies > 1) {
                            # make n_copies duplicates, so need n_copies - 1 added
                            current_label = labels[k]
                            to_replace = current_label
                            for (l in seq_len(n_copies - 1)) {
                                replacement = paste("(", current_label, "_", l, ":0,", current_label, "_", (l+1), ":0)", sep="")
                                str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                                to_replace = paste(current_label, "_", (l+1), sep="")
                            }
                        }
                        else {
                            # append _1 to label so that it matches new_suffices
                            to_replace = labels[k]
                            replacement = paste(to_replace, "_1", sep="")
                            str_newick_to_alter = gsub(to_replace, replacement, str_newick_to_alter)
                        }
                    }
                    new_obj@tumor_subclusters$hc[[sample_name]] = as.hclust(read.tree(text=str_newick_to_alter))
                    new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, (infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order]), drop=FALSE]
                    futur_colnames[(i:(i + n_cells - 1))] = paste(colnames(infercnv_obj@expr.data[, (infercnv_obj@observation_grouped_cell_indices[[sample_name]][new_data_order]), drop=FALSE]), new_suffixes, sep="_")
                }
                else {
                    if (length(infercnv_obj@tumor_subclusters$subclusters[[sample_name]]) > 1) {
                        flog.error(paste("No hclust object available for ", sample_name, " even though there is more than 1 cell", sep=""))
                        stop("Missing clustering information")
                    }

                    new_data_order = rep(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[1]], n_cells)
                    new_suffixes = seq_len(n_cells)

                    current_label = colnames(infercnv_obj@expr.data[, infercnv_obj@tumor_subclusters$subclusters[[sample_name]], drop=FALSE])
                    # str_newick_to_alter = paste(paste(lapply(seq_len((n_copies) - 1), function(l) {
                    #         paste("(", current_label, "_", l, ":0", sep="")
                    #     }), collapse=","), ",", current_label, "_", n_copies, ":0", paste(rep("):0", (n_copies - 2)), collapse=""), ");", sep="")
                    new_obj@tumor_subclusters$hc[[sample_name]] <- list()  # initialize empty object
                    new_obj@tumor_subclusters$hc[[sample_name]]$merge <- matrix(c(-1, -2), ncol=2, byrow=TRUE)
                    for (k in seq_len(n_copies - 2)) {
                        new_obj@tumor_subclusters$hc[[sample_name]]$merge <- rbind(new_obj@tumor_subclusters$hc[[sample_name]]$merge, c(k, -(k+2)))
                    }
                    new_obj@tumor_subclusters$hc[[sample_name]]$height <- rep(1, (n_copies - 1))    # define merge heights
                    new_obj@tumor_subclusters$hc[[sample_name]]$order <- seq_len(n_copies)              # order of leaves(trivial if hand-entered)
                    new_obj@tumor_subclusters$hc[[sample_name]]$labels <- paste(current_label, seq_len(n_copies), sep="_")    # lnew_obj@tumor_subclusters$hc[[sample_name]]bels of lenew_obj@tumor_subclusters$hc[[sample_name]]ves
                    class(new_obj@tumor_subclusters$hc[[sample_name]]) <- "hclust"
                    new_obj@expr.data[, (i:(i + n_cells - 1))] = infercnv_obj@expr.data[, new_data_order, drop=FALSE]
                    futur_colnames[(i:(i + n_cells - 1))] = paste(colnames(infercnv_obj@expr.data[, new_data_order, drop=FALSE]), new_suffixes, sep="_")
                }
            }
            # else {
                ## should never happen
            # }

            new_obj@tumor_subclusters$subclusters[[sample_name]][[paste(sample_name, "s1", sep="_")]] = c(i:(i + length(sampled_indices) - 1))
            new_obj@observation_grouped_cell_indices[[sample_name]]=c(i:(i + length(sampled_indices) - 1))
            i = i + length(sampled_indices)
        }
    }
    else { # do nothing on observations
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {
            new_obj@expr.data[, (i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))] = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE]
            new_obj@observation_grouped_cell_indices[[sample_name]] = (i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))
            new_obj@tumor_subclusters$hc[[sample_name]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]
            new_obj@tumor_subclusters$subclusters[[sample_name]] = infercnv_obj@tumor_subclusters$subclusters[[sample_name]]
            futur_colnames[(i:(i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]]) - 1))] = colnames(infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE])
            i = i + length(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
        }
    }
    colnames(new_obj@expr.data) = futur_colnames

    # now drop empty columns/cells that were made in case of every_n sampling missing some subclusters
    if (i <= dim(new_obj@expr.data)[2]) {
        new_obj@expr.data = new_obj@expr.data[, -(i:dim(new_obj@expr.data)[2]), drop=FALSE]
    }

    return(new_obj)
}


#### try to add support for k_obs_groups based splitting and not only cluster_by_groups=TRUE



#' @title plot_per_group
#' 
#' @description Takes an infercnv object and subdivides it into one object per group of cells 
#' to allow plotting of each group on a seperate plot. If references are selected, they will appear
#' on the observation heatmap area as it is larger.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param sample Whether unique groups of cells should be sampled from or not. (see other parameters for how sampling is done) (Default: FALSE)
#' 
#' @param n_cells Number of cells that should be sampled per group if sampling is enabled (default = 1000) .
#'
#' @param above_m Sample only groups that have at least above_m cells if sampling is enabled. (default: 1000)
#' Does not require every_n to be set.
#'
#' @param k_obs_groups Number of groups to break each group in with cutree (in the color bars on the left side of the plot only). (Default: 1)
#' 
#' @param every_n Sample 1 cell every_n cells for each group that has above_m cells, if sampling is enabled. 
#' If subclusters are defined, this will make sure that at least one cell per subcluster is sampled. 
#' Requires above_m to be set to work, overriding n_cells parameter. (Default: NULL)
#' 
#' @param base_filename Base prefix for the output files names. 
#' Will be followed by OBS/REF to indidate the type of the group, and the group name. (Default: "infercnv_per_group")
#'
#' @param output_format Output format for the figure. Choose between "png", "pdf" and NA. NA means to only write the text outputs without generating the figure itself. (default: "png")
#'
#' @param write_expr_matrix Includes writing a matrix file containing the expression data that is plotted in the heatmap. (default: FALSE)
#' 
#' @param save_objects Whether to save the infercnv objects generated for each group as RDS. (default: FALSE)
#'
#' @param png_res Resolution for png output. (Default: 300)
#'
#' @param dynamic_resize Factor (>= 0) by which to scale the dynamic resize of the observation 
#'                       heatmap and the overall plot based on how many cells there are.
#'                       Default is 0, which disables the scaling. Try 1 first if you want to enable. (Default: 0)
#'
#' @param out_dir Directory in which to save plots and other outputs.
#'
#' @param on_observations  boolean (default=TRUE), plot observations data (tumor cells).
#'
#' @param on_references  boolean (default=TRUE), plot references (normal cells).
#'
#' @param useRaster Whether to use rasterization for drawing heatmap. Only disable if it produces an error as it is much faster than not using it.
#' 
#' @return void
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
#' infercnv::plot_per_group(infercnv_object_example, out_dir=tempfile())
#'

plot_per_group <- function(infercnv_obj,
    on_references=TRUE,
    on_observations=TRUE,
    sample=FALSE,
    n_cells=1000,
    every_n=NULL,
    above_m=1000,
    k_obs_groups=1,
    base_filename="infercnv_per_group",
    output_format="png",
    write_expr_matrix=TRUE,
    save_objects=FALSE,
    png_res=300,
    dynamic_resize=0,
    useRaster=TRUE,
    out_dir) {

    plot_center = mean(infercnv_obj@expr.data)
    plot_range = quantile(infercnv_obj@expr.data[infercnv_obj@expr.data != plot_center], c(0.01, 0.99))

    if (on_references == TRUE) {
        for (sample_name in names(infercnv_obj@reference_grouped_cell_indices)) {

            new_obj <- new(
                Class = "infercnv",
                expr.data = infercnv_obj@expr.data[, infercnv_obj@reference_grouped_cell_indices[[sample_name]], drop=FALSE], 
                count.data = matrix(),  # not needed
                gene_order = infercnv_obj@gene_order,
                # reference_grouped_cell_indices = list(sample_name = seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]])),
                # observation_grouped_cell_indices = list(),
                reference_grouped_cell_indices = list(),
                observation_grouped_cell_indices = list(),
                tumor_subclusters = list(),
                .hspike = infercnv_obj@.hspike)

            new_obj@observation_grouped_cell_indices[[sample_name]] = seq_along(infercnv_obj@reference_grouped_cell_indices[[sample_name]])

            if (!is.null(infercnv_obj@tumor_subclusters)) {
                new_obj@tumor_subclusters$hc = list()
                new_obj@tumor_subclusters$subclusters = list()
                new_obj@tumor_subclusters$hc[["all_observations"]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]

                # match(pre_obj@observation_grouped_cell_indices$CBTP3_187188XL, unlist(pre_obj@tumor_subclusters$subclusters$CBTP3_187188XL, use.names=FALSE))
                # match -> pos1 in pos2
                for (subcluster_id in names(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])) {
                    new_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster_id]] = unlist(match(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]], infercnv_obj@reference_grouped_cell_indices[[sample_name]]))
                }
            }
            else {
                new_obj@tumor_subclusters = NULL
            }

            if (sample) {
                if (!is.null(above_m) && ncol(new_obj@expr.data) > above_m) {
                    new_obj <- sample_object(new_obj,
                                             n_cells=n_cells,
                                             every_n=every_n,
                                             above_m=above_m,
                                             on_references=FALSE,
                                             on_observations=TRUE)
                }
            }

            if (save_objects) {
                saveRDS(new_obj, file.path(out_dir, paste(base_filename, "_REF_", make_filename(sample_name), ".infercnv_obj", sep="")))
            }

            plot_cnv(new_obj,
                out_dir=out_dir,
                title=paste("inferCNV", sample_name),
                obs_title=sample_name,
                ref_title="",
                cluster_by_groups=FALSE,
                cluster_references=FALSE,
                k_obs_groups=k_obs_groups,
                contig_cex=1,
                x.center=plot_center,
                x.range=plot_range,
                color_safe_pal=FALSE,
                output_filename=paste(base_filename, "REF", make_filename(sample_name), sep="_"),
                output_format=output_format, #pdf, png, NA
                png_res=png_res,
                dynamic_resize=dynamic_resize,
                ref_contig=NULL,
                write_expr_matrix=write_expr_matrix,
                useRaster=useRaster
                )
        }
    }

    if (on_observations == TRUE) {
        for (sample_name in names(infercnv_obj@observation_grouped_cell_indices)) {

            new_obj <- new(
                Class = "infercnv",
                expr.data = infercnv_obj@expr.data[, infercnv_obj@observation_grouped_cell_indices[[sample_name]], drop=FALSE], 
                count.data = matrix(),  # not needed
                gene_order = infercnv_obj@gene_order,
                reference_grouped_cell_indices = list(),
                observation_grouped_cell_indices = list(),
                tumor_subclusters = list(),
                .hspike = infercnv_obj@.hspike)

            new_obj@observation_grouped_cell_indices[[sample_name]] = seq_along(infercnv_obj@observation_grouped_cell_indices[[sample_name]])
            if (!is.null(infercnv_obj@tumor_subclusters)) {
                new_obj@tumor_subclusters$hc = list()
                new_obj@tumor_subclusters$subclusters = list()
                new_obj@tumor_subclusters$hc[["all_observations"]] = infercnv_obj@tumor_subclusters$hc[[sample_name]]
             
                for (subcluster_id in names(infercnv_obj@tumor_subclusters$subclusters[[sample_name]])) {
                    new_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster_id]] = unlist(match(infercnv_obj@tumor_subclusters$subclusters[[sample_name]][[subcluster_id]], infercnv_obj@observation_grouped_cell_indices[[sample_name]]))
                    names(new_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster_id]]) = colnames(new_obj@expr.data[, new_obj@tumor_subclusters$subclusters[["all_observations"]][[subcluster_id]], drop=FALSE])
                }
            }
            else {
                new_obj@tumor_subclusters = NULL
            }


            if (sample) {
                if (!is.null(above_m) && ncol(new_obj@expr.data) > above_m) {
                    new_obj <- sample_object(new_obj,
                                             n_cells=n_cells,
                                             every_n=every_n,
                                             above_m=above_m,
                                             on_references=FALSE,
                                             on_observations=TRUE)
                }
            }

            if (save_objects) {
                saveRDS(new_obj, file.path(out_dir, paste(base_filename, "_OBS_", make_filename(sample_name), ".infercnv_obj", sep="")))
            }

            plot_cnv(new_obj,
                out_dir=out_dir,
                title=paste("inferCNV", sample_name),
                obs_title=sample_name,
                ref_title="",
                cluster_by_groups=FALSE,
                cluster_references=FALSE,
                k_obs_groups=k_obs_groups,
                contig_cex=1,
                x.center=plot_center,
                x.range=plot_range,
                color_safe_pal=FALSE,
                output_filename=paste(base_filename, "OBS", make_filename(sample_name), sep="_"),
                output_format=output_format, #pdf, png, NA
                png_res=png_res,
                dynamic_resize=dynamic_resize,
                ref_contig=NULL,
                write_expr_matrix=write_expr_matrix,
                useRaster=useRaster
                )
        }
    }
}

make_filename <- function(text) {
    text <- gsub("[/\\:*?\"<>|]", "_", text)
}

