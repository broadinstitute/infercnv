
#' Function doing the actual analysis before calling the plotting functions.
#'
#' @title run() : Invokes a routine inferCNV analysis to Infer CNV changes given a matrix of RNASeq counts.
#'
#' @param infercnv_obj An infercnv object populated with raw count data
#'
#' @param cutoff Cut-off for the min average read counts per gene among reference cells. (default: 1)
#'
#' @param min_cells_per_gene minimum number of reference cells requiring expression measurements to include the corresponding gene.
#'                           default: 3
#'
#' @param out_dir path to directory to deposit outputs (default: NULL, required to provide non NULL)
#'
#' ## Smoothing params
#' @param window_length Length of the window for the moving average
#'                          (smoothing). Should be an odd integer. (default: 101)#'
#'
#' @param smooth_method  Method to use for smoothing: c(runmeans,pyramidinal,coordinates)  default: pyramidinal
#'
#' #####
#'
#' @param num_ref_groups The number of reference groups or a list of
#'                           indices for each group of reference indices in
#'                           relation to reference_obs. (default: NULL)
#'
#' @param ref_subtract_use_mean_bounds   Determine means separately for each ref group, then remove intensities within bounds of means (default: TRUE)
#'                                       Otherwise, uses mean of the means across groups.
#'
#' #############################
#'
#' @param cluster_by_groups   If observations are defined according to groups (ie. patients), each group
#'                            of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
#'
#' @param cluster_references Whether to cluster references within their annotations or not. (dendrogram not displayed)
#'                             (default: TRUE)
#'
#'
#' @param k_obs_groups Number of groups in which to break the observations. (default: 1)
#'
#'
#'
#' @param hclust_method Method used for hierarchical clustering of cells. Valid choices are:
#' "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid".
#' default("ward.D2")
#'
#' @param max_centered_threshold The maximum value a value can have after
#'                                   centering. Also sets a lower bound of
#'                                   -1 * this value. (default: 3),
#'                               can set to a numeric value or "auto" to bound by the mean bounds across cells.
#'                               Set to NA to turn off.
#'
#' @param scale_data  perform Z-scaling of logtransformed data (default: FALSE).  This may be turned on if you have
#'                    very different kinds of data for your normal and tumor samples. For example, you need to use GTEx
#'                    representative normal expression profiles rather than being able to leverage normal single cell data
#'                    that goes with your experiment.
#'
#' #########################################################################
#' ## Downstream Analyses (HMM or non-DE-masking) based on tumor subclusters
#'
#' @param HMM  when set to True, runs HMM to predict CNV level (default: FALSE)
#'
#' @param HMM_transition_prob   transition probability in HMM (default: 1e-6)
#'
#'
#' @param HMM_report_by   cell, consensus, subcluster (default: subcluster)  Note, reporting is performed entirely separately from the HMM prediction.  So, you can predict on subclusters, but get per-cell level reporting (more voluminous output).
#'
#'
#' @param HMM_type  HMM model type. Options: (i6 or i3):
#'                          i6: infercnv 6-state model (0, 0.5, 1, 1.5, 2, >2) where state emissions are calibrated based on simulated CNV levels.
#'                          i3: infercnv 3-state model (del, neutral, amp) configured based on normal cells and HMM_i3_pval
#'
#' @param HMM_i3_pval     p-value for HMM i3 state overlap (default: 0.05)
#'
#'
#' @param HMM_i3_use_KS   boolean: use the KS test statistic to estimate mean of amp/del distributions (ala HoneyBadger). (default=TRUE)
#'
#'
#' ## Filtering low-conf HMM preds via BayesNet P(Normal)
#'
#' @param BayesMaxPNormal  maximum P(Normal) allowed for a CNV prediction according to BayesNet. (default=0.5, note zero turns it off)
#' 
#' @param reassignCNVs (boolean) Given the CNV associated probability of belonging to each possible state, 
#'                      reassign the state assignments made by the HMM to the state that has the highest probability. (default: TRUE)
#'
#' ######################
#' ## Tumor subclustering
#'
#' @param analysis_mode options(samples|subclusters|cells), Grouping level for image filtering or HMM predictions.
#'                                                          default: samples (fastest, but subclusters is ideal)
#'
#' @param tumor_subcluster_partition_method  method for defining tumor subclusters. Options('leiden', 'random_trees', 'qnorm')
#'                                           leiden: Runs a nearest neighbor search, where communities are then partitionned with the Leiden algorithm.
#'                                           random_trees: Slow, uses permutation statistics w/ tree construction.
#'                                           qnorm: defines tree height based on the quantile defined by the tumor_subcluster_pval
#'
#' @param tumor_subcluster_pval max p-value for defining a significant tumor subcluster (default: 0.1)
#'
#' @param k_nn number k of nearest neighbors to search for when using the Leiden partition method for subclustering (default: 20)
#'
#' @param leiden_method Method used to generate the graph on which the Leiden algorithm is applied, one of "PCA" or "simple". (default: "PCA")
#'
#' @param leiden_function Whether to use the Constant Potts Model (CPM) or modularity in igraph. Must be either "CPM" or "modularity". (default: "CPM")
#'
#' @param leiden_resolution resolution parameter for the Leiden algorithm using the CPM quality score (default: auto)
#'
#' @param leiden_method_per_chr Method used to generate the graph on which the Leiden algorithm is applied for the per chromosome subclustering, one of "PCA" or "simple". (default: "simple")
#'
#' @param leiden_function_per_chr Whether to use the Constant Potts Model (CPM) or modularity in igraph for the per chromosome subclustering. Must be either "CPM" or "modularity". (default: "modularity")
#'
#' @param leiden_resolution_per_chr resolution parameter for the Leiden algorithm for the per chromosome subclustering (default: 1)
#'
#' @param per_chr_hmm_subclusters Run subclustering per chromosome over all cells combined to run the HMM on those subclusters instead. Only applicable when using Leiden subclustering.
#'                                This should provide enough definition in the predictions while avoiding subclusters that are too small thus providing less evidence to work with. (default: FALSE)
#'
#' @param per_chr_hmm_subclusters_references Whether the per chromosome subclustering should also be done on references, which should not have as much variation as observations. (default = FALSE)
#'
#' 
#' @param z_score_filter Z-score used as a treshold to filter genes used for subclustering. 
#'                       Applied based on reference genes to automatically ignore genes with high expression variability such as MHC genes. (default: 0.8)
#'
#'
#' #############################
#' ## de-noising parameters ####
#'
#' @param denoise       If True, turns on denoising according to options below
#'
#' @param noise_filter  Values +- from the reference cell mean will be set to zero (whitening effect)
#'                      default(NA, instead will use sd_amplifier below.
#'
#' @param sd_amplifier  Noise is defined as mean(reference_cells) +- sdev(reference_cells) * sd_amplifier
#'                      default: 1.5
#'
#' @param noise_logistic use the noise_filter or sd_amplifier based threshold (whichever is invoked) as the midpoint in a
#'                       logistic model for downscaling values close to the mean. (default: FALSE)
#'
#'
#' ##################
#' ## Outlier pruning
#'
#' @param outlier_method_bound Method to use for bounding outlier values. (default: "average_bound")
#'                             Will preferentially use outlier_lower_bounda and outlier_upper_bound if set.
#' @param outlier_lower_bound  Outliers below this lower bound will be set to this value.
#' @param outlier_upper_bound  Outliers above this upper bound will be set to this value.
#'
#'
#' ##########################
#' ## Misc options
#'
#' @param final_scale_limits The scale limits for the final heatmap output by the run() method. Default "auto". Alt, c(low,high)
#'
#' @param final_center_val   Center value for final heatmap output by the run() method.
#'
#' @param debug If true, output debug level logging.
#'
#' @param num_threads (int) number of threads for parallel steps (default: 4)
#'
#' @param plot_steps If true, saves infercnv objects and plots data at the intermediate steps.
#'
#' @param inspect_subclusters If true, plot subclusters as annotations after the subclustering step to easily see if the subclustering options are good. (default = TRUE)
#'
#' @param resume_mode  leverage pre-computed and stored infercnv objects where possible. (default=TRUE)
#'
#' @param png_res Resolution for png output.
#'
#' @param no_plot   don't make any of the images. Instead, generate all non-image outputs as part of the run. (default: FALSE)
#'
#' @param no_prelim_plot  don't make the preliminary infercnv image (default: FALSE)
#'
#' @param write_expr_matrix Whether to write text files with the content of matrices when generating plots (default: FALSE)
#'
#' @param write_phylo Whether to write newick strings of the dendrograms displayed on the left side of the heatmap to file (default: FALSE)
#'
#' @param output_format Output format for the figure. Choose between "png", "pdf" and NA. NA means to only write the text outputs without generating the figure itself. (default: "png")
#'
#' @param plot_chr_scale Whether to scale the chromosme width on the heatmap based on their actual size rather than just the number of expressed genes.
#'
#' @param chr_lengths A named list of chromsomes lengths to use when plot_chr_scale=TRUE, or else chromosome size is assumed to be the last chromosome's stop position + 10k bp
#'
#' @param useRaster Whether to use rasterization for drawing heatmap. Only disable if it produces an error as it is much faster than not using it. (default: TRUE)
#'
#' @param plot_probabilities option to plot posterior probabilities (default: TRUE)
#'
#' @param save_rds Whether to save the current step object results as an .rds file (default: TRUE)
#'
#' @param save_final_rds Whether to save the final object results as an .rds file (default: TRUE)
#'
#' @param diagnostics option to create diagnostic plots after running the Bayesian model (default: FALSE)
#'
#' #######################
#' ## Experimental options
#'
#' @param remove_genes_at_chr_ends experimental option: If true, removes the window_length/2 genes at both ends of the chromosome.
#'
#' @param prune_outliers  Define outliers loosely as those that exceed the mean boundaries among all cells.  These are set to the bounds.
#'
#' ## experimental opts involving DE analysis
#'
#' @param mask_nonDE_genes If true, sets genes not significantly differentially expressed between tumor/normal to
#'                          the mean value for the complete data set (default: 0.05)
#'
#' @param mask_nonDE_pval  p-value threshold for defining statistically significant DE genes between tumor/normal
#
#' @param test.use statistical test to use.  (default: "wilcoxon") alternatives include 'perm' or 't'.'
#'
#' @param require_DE_all_normals If mask_nonDE_genes is set, those genes will be masked only if they are are found as DE according to test.use and mask_nonDE_pval in each of the comparisons to normal cells options: {"any", "most", "all"} (default: "any")
#'
#' other experimental opts
#'
#' @param sim_method   method for calibrating CNV levels in the i6 HMM (default: 'meanvar')
#'
#' @param sim_foreground  don't use... for debugging, developer option.
#'
#' @param hspike_aggregate_normals  instead of trying to model the different normal groupings individually, just merge them in the hspike.
#'
#' @param up_to_step run() only up to this exact step number (default: 100 >> 23 steps currently in the process)
#'
#' @return infercnv_obj containing filtered and transformed data
#'
#' @export
#'
#' @examples
#' data(infercnv_data_example)
#' data(infercnv_annots_example)
#' data(infercnv_genes_example)
#'
#' infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_example, 
#'                                                           gene_order_file=infercnv_genes_example,
#'                                                           annotations_file=infercnv_annots_example,
#'                                                           ref_group_names=c("normal"))
#'
#' infercnv_object_example <- infercnv::run(infercnv_object_example,
#'                                          cutoff=1,
#'                                          out_dir=tempfile(), 
#'                                          cluster_by_groups=TRUE, 
#'                                          denoise=TRUE,
#'                                          HMM=FALSE,
#'                                          num_threads=2,
#'                                          analysis_mode="samples",
#'                                          no_plot=TRUE)
#'
run <- function(infercnv_obj,

                # gene filtering settings
                cutoff=1,
                min_cells_per_gene=3,

                out_dir=NULL,

                ## smoothing params
                window_length=101,
                smooth_method=c('pyramidinal', 'runmeans', 'coordinates'),

                num_ref_groups=NULL,
                ref_subtract_use_mean_bounds=TRUE,

                # observation cell clustering settings
                cluster_by_groups=TRUE,
                cluster_references=TRUE,
                k_obs_groups=1,

                hclust_method='ward.D2',

                max_centered_threshold=3, # or set to a specific value or "auto", or NA to turn off
                scale_data=FALSE,

                ## HMM opts
                HMM=FALSE, # turn on to auto-run the HMM prediction of CNV levels
                ## tumor subclustering opts
                HMM_transition_prob=1e-6,
                HMM_report_by=c("subcluster","consensus","cell"),
                HMM_type=c('i6', 'i3'),
                HMM_i3_pval=0.05,
                HMM_i3_use_KS=FALSE,
                BayesMaxPNormal=0.5,

                ## some experimental params
                #sim_method=c('meanvar', 'simple', 'splatter'), ## only meanvar supported, others experimental
                sim_method='meanvar',
                sim_foreground=FALSE, ## experimental
                reassignCNVs=TRUE,


                ## tumor subclustering options
                analysis_mode=c('subclusters', 'samples', 'cells'), # for filtering and HMM
                tumor_subcluster_partition_method=c('leiden', 'random_trees', 'qnorm', 'pheight', 'qgamma', 'shc'),
                tumor_subcluster_pval=0.1,
                k_nn=20,
                leiden_method=c("PCA", "simple"),
                leiden_function = c("CPM", "modularity"),
                leiden_resolution="auto",
                leiden_method_per_chr=c("simple", "PCA"),
                leiden_function_per_chr = c("modularity", "CPM"),
                leiden_resolution_per_chr = 1,
                per_chr_hmm_subclusters=FALSE,
                per_chr_hmm_subclusters_references=FALSE,
                z_score_filter = 0.8,


                ## noise settings
                denoise=FALSE,
                noise_filter=NA,
                sd_amplifier = 1.5,
                noise_logistic=FALSE, # if false, does complete 'noise' elimination.

                # outlier adjustment settings
                outlier_method_bound="average_bound",
                outlier_lower_bound=NA,
                outlier_upper_bound=NA,

                ## misc options
                final_scale_limits = NULL,
                final_center_val = NULL,
                debug=FALSE, #for debug level logging
                num_threads = 4,
                plot_steps=FALSE,
                inspect_subclusters = TRUE,
                resume_mode=TRUE,
                png_res=300,
                plot_probabilities = TRUE,
                save_rds = TRUE,
                save_final_rds = TRUE,
                diagnostics = FALSE,

                ## experimental options
                remove_genes_at_chr_ends=FALSE,
                prune_outliers=FALSE,

                mask_nonDE_genes=FALSE,
                mask_nonDE_pval=0.05, # use permissive threshold
                test.use='wilcoxon',
                require_DE_all_normals="any",


                hspike_aggregate_normals = FALSE,

                no_plot = FALSE,
                no_prelim_plot = FALSE,
                write_expr_matrix = FALSE,
                write_phylo = FALSE,
                output_format = "png",
                plot_chr_scale = FALSE,
                chr_lengths = NULL,
                useRaster = TRUE,

                up_to_step=100

) {


    smooth_method = match.arg(smooth_method)
    HMM_type = match.arg(HMM_type)
    if (HMM && HMM_type == 'i6' && smooth_method == 'coordinates') {
        flog.error("Cannot use 'coordinate' smoothing method with i6 HMM model at this time.")
        stop("Incompatible HMM mode and smoothing method.")
    }
    if (smooth_method == 'coordinates' && window_length < 10000) {
        window_length=10000000
        flog.warn(paste0("smooth_method set to 'coordinates', but window_length ",
                         "is less than 10.000, setting it to 10.000.000. Please ",
                         "set a different value > 10.000 if this default does not seem suitable."))
    }

    leiden_function = match.arg(leiden_function)
    leiden_function_per_chr = match.arg(leiden_function_per_chr)
    HMM_report_by = match.arg(HMM_report_by)
    analysis_mode = match.arg(analysis_mode)
    if(analysis_mode == "cells" && HMM_report_by != "cell") {
        flog.warn(paste0("analysis_mode is \"cells\" but HMM_report_by is not, "),
                         "changing HMM_report_by to \"cells\".")
        HMM_report_by = "cell"
    }
    tumor_subcluster_partition_method = match.arg(tumor_subcluster_partition_method)
    leiden_method = match.arg(leiden_method)
    leiden_method_per_chr = match.arg(leiden_method_per_chr)
    if (tumor_subcluster_partition_method != "leiden") {
        per_chr_hmm_subclusters = FALSE
    }
    
    if (debug) {
        flog.threshold(DEBUG)
    } else {
        flog.threshold(INFO)
    }
    
    flog.info(paste("::process_data:Start", sep=""))
    
    infercnv.env$GLOBAL_NUM_THREADS <- num_threads
    if (is.null(out_dir)) {
        flog.error("Error, out_dir is NULL, please provide a path.")
        stop("out_dir is NULL")
    }

    call_match = match.call()

    arg_names = names(call_match)
    arg_names = arg_names[-which(arg_names == "" | arg_names == "infercnv_obj")]

    # evaluate variables such as output_dir
    for (n in arg_names) {
        call_match[[n]] = get(n)
    }

    # add the arguments that are checked against a list of allowed values and autofilled with the first if none is given using match.arg()
    call_match$smooth_method = smooth_method
    call_match$HMM_type = HMM_type
    call_match$HMM_report_by = HMM_report_by
    call_match$analysis_mode = analysis_mode
    call_match$tumor_subcluster_partition_method = tumor_subcluster_partition_method
    call_match$hclust_method = hclust_method
    call_match$leiden_function = leiden_function
    # add argument needed to get relevant args list
    call_match$HMM = HMM
    call_match$z_score_filter = z_score_filter
    call_match$tumor_subcluster_pval = tumor_subcluster_pval
    call_match$k_nn = k_nn
    call_match$leiden_method = leiden_method
    call_match$per_chr_hmm_subclusters = per_chr_hmm_subclusters
    call_match$cluster_by_groups = cluster_by_groups
    call_match$max_centered_threshold = max_centered_threshold
    call_match$remove_genes_at_chr_ends = remove_genes_at_chr_ends
    call_match$prune_outliers = prune_outliers
    call_match$BayesMaxPNormal = BayesMaxPNormal
    call_match$mask_nonDE_genes = mask_nonDE_genes
    call_match$denoise = denoise
    call_match$noise_filter = noise_filter
    call_match$sd_amplifier = sd_amplifier
    call_match$noise_logistic = noise_logistic

    call_match$infercnv_obj = NULL

    # pull the list of argument that are set by the user in case they differ from the default
    current_args = as.list(call_match[-which(names(call_match) == "" | names(call_match) == "infercnv_obj")])
    # add the arguments that are checked against a list of allowed values and autofilled with the first if none is given using match.arg()
    # current_args$smooth_method = smooth_method
    # current_args$HMM_type = HMM_type
    # current_args$HMM_report_by = HMM_report_by
    # current_args$analysis_mode = analysis_mode
    # current_args$tumor_subcluster_partition_method = tumor_subcluster_partition_method

    if(out_dir != "." && !file.exists(out_dir)){
        flog.info(paste0("Creating output path ", out_dir))
        dir.create(out_dir)
    }

    infercnv_obj@options = c(current_args, infercnv_obj@options)

    #run_call <- match.call()
    call_match[[1]] <- as.symbol(".get_relevant_args_list")
    reload_info = eval(call_match)

    # reload_info$relevant_args
    # reload_info$expected_file_names

    skip_hmm = 0
    # skip_mcmc = 0
    skip_past = 0
    
    if (resume_mode) {
        flog.info("Checking for saved results.")
        for (i in rev(seq_along(reload_info$expected_file_names))) {
            if (file.exists(reload_info$expected_file_names[[i]])) {
                flog.info(paste0("Trying to reload from step ", i))
                # if ((i == 17 && skip_hmm == 0) || (i %in% (18:20))) {
                if ((i == 20) || (i %in% seq(17, 19) && skip_hmm == 0) || (i == 17 && skip_hmm == 2)) {   # step 17 appears in two conditions because if last found step is 18, step 19 requires the results of both step 17 and 18 to run
                    if (i == 18 && BayesMaxPNormal > 0) {  # mcmc_obj
                        mcmc_obj = readRDS(reload_info$expected_file_names[[i]])
                        if (!.compare_args(infercnv_obj@options, unlist(reload_info$relevant_args[1:i]), mcmc_obj@options)) {
                            rm(mcmc_obj)
                            invisible(gc())
                        }
                        else {
                            mcmc_obj@options = infercnv_obj@options
                            skip_hmm = max(skip_hmm, i - 16)
                            flog.info(paste0("Using backup MCMC from step ", i))
                        }
                    }
                    else {
                        hmm.infercnv_obj = readRDS(reload_info$expected_file_names[[i]])
                        if (!.compare_args(infercnv_obj@options, unlist(reload_info$relevant_args[1:i]), hmm.infercnv_obj@options)) {
                            rm(hmm.infercnv_obj)
                            invisible(gc())
                        }
                        else {
                            hmm.infercnv_obj@options = infercnv_obj@options
                            skip_hmm = max(skip_hmm, i - 16)  # max for case where step 18 is found, so can be skipped, but still needed to reload step 17 results to be able to run step 19
                            flog.info(paste0("Using backup HMM from step ", i))
                        }
                    }
                }
                else if (!(i %in% 17:20)) {
                    reloaded_infercnv_obj = readRDS(reload_info$expected_file_names[[i]])
                    if (skip_past > i) { # in case denoise was found
                        if (21 > i) { # if 22/21 already found and checked HMM too, stop
                            break
                        }
                        else { # if 22 denoise found, don't check 21 maskDE
                            next
                        }
                    }
                    if (.compare_args(infercnv_obj@options, unlist(reload_info$relevant_args[1:i]), reloaded_infercnv_obj@options)) {
                        if (i ==  15 && per_chr_hmm_subclusters) {
                            tmp = nchar(reload_info$expected_file_names[[i]])
                            tmp = paste0(substr(reload_info$expected_file_names[[i]], 1, tmp-13), ".per_chr_subclusters", substr(reload_info$expected_file_names[[i]], tmp-12, tmp))
                            if (file.exists(tmp)) {
                                subclusters_per_chr = readRDS(tmp)
                            }
                            else {
                                next
                            }
                        }
                        options_backup = infercnv_obj@options
                        infercnv_obj = reloaded_infercnv_obj # replace input infercnv_obj
                        rm(reloaded_infercnv_obj) # remove first (temporary) reference so there's no duplication when they would diverge
                        infercnv_obj@options = options_backup
                        skip_past = i
                        flog.info(paste0("Using backup from step ", i))
                        if (i < 21) { # do not stop check right away if denoise was found to also allow to check for HMM
                            break
                        }
                    }
                    else {
                        rm(reloaded_infercnv_obj)
                        invisible(gc())
                    }
                }
            }
        }
    }
 

    step_count = 0;
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 1
    flog.info(sprintf("\n\n\tSTEP %d: incoming data\n", step_count))

    if (skip_past < step_count && save_rds) {
        saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
    }
    
    ## #################################################
    ## Step: removing insufficiently expressed genes
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 2
    flog.info(sprintf("\n\n\tSTEP %02d: Removing lowly expressed genes\n", step_count))
    
    ## Remove genes that aren't sufficiently expressed, according to min mean count cutoff.
    ## Examines the original (non-log-transformed) data, gets mean for each gene, and removes genes
    ##  with mean values below cutoff.

    if (skip_past < step_count) {
    # }
        infercnv_obj <- require_above_min_mean_expr_cutoff(infercnv_obj, cutoff)
        
        ## require each gene to be present in a min number of cells for ref sets
        
        infercnv_obj <- require_above_min_cells_ref(infercnv_obj, min_cells_per_gene=min_cells_per_gene)
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
    }
    
    
    ## #########################################
    ## # STEP: normalization by sequencing depth

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 3
    flog.info(sprintf("\n\n\tSTEP %02d: normalization by sequencing depth\n", step_count))
    
    
    resume_file_token = ifelse( (HMM), paste0("HMM",HMM_type), "")

    if (skip_past < step_count) {
        infercnv_obj <- normalize_counts_by_seq_depth(infercnv_obj)
        
        if (HMM && HMM_type == 'i6') {
            ## add in the hidden spike needed by the HMM
            infercnv_obj <- .build_and_add_hspike(infercnv_obj, sim_method=sim_method, aggregate_normals=hspike_aggregate_normals)
            
            if (sim_foreground) {
                infercnv_obj <- .sim_foreground(infercnv_obj, sim_method=sim_method)
            }
        }
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
    }
    

    ## #########################
    ## Step: log transformation

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 4
    flog.info(sprintf("\n\n\tSTEP %02d: log transformation of data\n", step_count))

    if (skip_past < step_count) {
        infercnv_obj <- log2xplus1(infercnv_obj)

        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())

        ## Plot incremental steps.
        if (plot_steps){
            plot_cnv(infercnv_obj=infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_log_transformed_data",step_count),
                     output_filename=sprintf("infercnv.%02d_log_transformed",step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster
                     )
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 5
    if (scale_data) {
        
        flog.info(sprintf("\n\n\tSTEP %02d: scaling all expression data\n", step_count))

        if (skip_past < step_count) {
            infercnv_obj <- scale_infercnv_expr(infercnv_obj)
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_scaled",step_count),
                         output_filename=sprintf("infercnv.%02d_scaled",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    ## #################################################
    ## Step: Split the reference data into groups if requested
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 6
    if (!is.null(num_ref_groups)) {
        
        if (! has_reference_cells(infercnv_obj)) {
            stop("Error, no reference cells defined. Cannot split them into groups as requested")
        }
        
        flog.info(sprintf("\n\n\tSTEP %02d: splitting reference data into %d clusters\n", step_count, num_ref_groups))

        if (skip_past < step_count) {
            infercnv_obj <- split_references(infercnv_obj,
                                       num_groups=num_ref_groups,
                                       hclust_method=hclust_method)
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
        }
        
    }
    
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 7
    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method == 'random_trees') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
        
        resume_file_token = paste0(resume_file_token, ".rand_trees")

        if (skip_past < step_count) {
            infercnv_obj <- define_signif_tumor_subclusters_via_random_smooothed_trees(infercnv_obj,
                                                                                       p_val=tumor_subcluster_pval,
                                                                                       hclust_method=hclust_method,
                                                                                       cluster_by_groups=cluster_by_groups)
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (plot_steps) {
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         output_filename=sprintf("infercnv.%02d_tumor_subclusters.%s", step_count, tumor_subcluster_partition_method),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }

            if (inspect_subclusters) {
                plot_subclusters(infercnv_obj,
                                 out_dir=out_dir,
                                 output_filename="infercnv_subclusters")
            }
        }
    }
    
    
    ## ##################################
    ## Step: Subtract average reference
    ## Since we're in log space, this now becomes log(fold_change)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 8
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (before smoothing)\n", step_count))

    if (skip_past < step_count) {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("infercnv.%02d_remove_average", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 9
    if (! is.na(max_centered_threshold)) {
        
        ## #####################################################
        ## Apply maximum centered expression thresholds to data
        ## Cap values between threshold and -threshold, retaining earlier center_cell_expr_across_chromosome
        
        flog.info(sprintf("\n\n\tSTEP %02d: apply max centered expression threshold: %s\n", step_count, max_centered_threshold))

        if (skip_past < step_count) {
            threshold = max_centered_threshold
            if (is.character(max_centered_threshold) && max_centered_threshold == "auto") {
                threshold = mean(abs(get_average_bounds(infercnv_obj)))
                flog.info(sprintf("Setting max centered threshoolds via auto to: +- %g", threshold))
            }
            
            infercnv_obj <- apply_max_threshold_bounds(infercnv_obj, threshold=threshold)
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_apply_max_centered_expr_threshold",step_count),
                         output_filename=sprintf("infercnv.%02d_apply_max_centred_expr_threshold",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    
    ## #########################################################################
    ## Step: For each cell, smooth the data along chromosome with gene windows
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 10
    flog.info(sprintf("\n\n\tSTEP %02d: Smoothing data per cell by chromosome\n", step_count))
    
    if (skip_past < step_count) {

        if (smooth_method == 'runmeans') {
            
            infercnv_obj <- smooth_by_chromosome_runmeans(infercnv_obj, window_length)
        } else if (smooth_method == 'pyramidinal') {
            
            infercnv_obj <- smooth_by_chromosome(infercnv_obj, window_length=window_length, smooth_ends=TRUE)
        } else if (smooth_method == 'coordinates') {
            infercnv_obj <- smooth_by_chromosome_coordinates(infercnv_obj, window_length=window_length)
        } else {
            stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
        }
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        ## Plot incremental steps.
        if (plot_steps){
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_smoothed_by_chr",step_count),
                     output_filename=sprintf("infercnv.%02d_smoothed_by_chr", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    ##
    ## Step:
    ## Center cells/observations after smoothing. This helps reduce the
    ## effect of complexity.
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 11
    flog.info(sprintf("\n\n\tSTEP %02d: re-centering data across chromosome after smoothing\n", step_count))
    
    if (skip_past < step_count) {
        infercnv_obj <- center_cell_expr_across_chromosome(infercnv_obj, method="median")
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        ## Plot incremental steps.
        if (plot_steps) {
            
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_centering_of_smoothed",step_count),
                     output_filename=sprintf("infercnv.%02d_centering_of_smoothed", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            
        }
    }
    
    
    
    ## ##################################
    ## Step: Subtract average reference (adjustment)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 12
    flog.info(sprintf("\n\n\tSTEP %02d: removing average of reference data (after smoothing)\n", step_count))
    
    if (skip_past < step_count) {
        infercnv_obj <- subtract_ref_expr_from_obs(infercnv_obj, inv_log=FALSE, use_bounds=ref_subtract_use_mean_bounds)
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_remove_average",step_count),
                     output_filename=sprintf("infercnv.%02d_remove_average", step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
        }
    }
    
    
    ## Step: Remove Ends
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 13
    if (remove_genes_at_chr_ends == TRUE && smooth_method != 'coordinates') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: removing genes at chr ends\n", step_count))
        
        if (skip_past < step_count) {
            infercnv_obj <- remove_genes_at_ends_of_chromosomes(infercnv_obj, window_length)
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps){
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_remove_genes_at_chr_ends",step_count),
                         output_filename=sprintf("infercnv.%02d_remove_genes_at_chr_ends",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res)
                
            }
        }
    }
    
    
    ## ###########################
    ## Step: invert log transform  (convert from log(FC) to FC)
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 14
    flog.info(sprintf("\n\n\tSTEP %02d: invert log2(FC) to FC\n", step_count))

    if (skip_past < step_count) {
        
        infercnv_obj <- invert_log2(infercnv_obj)
        
        if (save_rds) {
            saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
        }
        invisible(gc())
        
        if (plot_steps) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title=sprintf("%02d_invert_log_transform log(FC)->FC",step_count),
                     output_filename=sprintf("infercnv.%02d_invert_log_FC",step_count),
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            
        }
    }
    
    
    ## ###################################################################
    ## Done restoring infercnv_obj's from files now under resume_mode
    ## ###################################################################

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 15
    if (analysis_mode == 'subclusters' & tumor_subcluster_partition_method != 'random_trees') {
        
        resume_file_token = paste0(resume_file_token, '.', tumor_subcluster_partition_method)

        flog.info(sprintf("\n\n\tSTEP %02d: computing tumor subclusters via %s\n", step_count, tumor_subcluster_partition_method))
        
        if (skip_past < step_count) {
            res <- define_signif_tumor_subclusters(infercnv_obj = infercnv_obj,
                                                   p_val = tumor_subcluster_pval,
                                                   k_nn = k_nn,
                                                   leiden_resolution = leiden_resolution,
                                                   leiden_method = leiden_method,
                                                   leiden_function = leiden_function,
                                                   leiden_method_per_chr = leiden_method_per_chr,
                                                   leiden_function_per_chr = leiden_function_per_chr,
                                                   leiden_resolution_per_chr = leiden_resolution_per_chr,
                                                   hclust_method = hclust_method,
                                                   cluster_by_groups = cluster_by_groups,
                                                   partition_method = tumor_subcluster_partition_method,
                                                   per_chr_hmm_subclusters = per_chr_hmm_subclusters,
                                                   z_score_filter=z_score_filter
                                                   )

            infercnv_obj = res[[1]]
            subclusters_per_chr = res[[2]]
            rm(res)

            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
                if(per_chr_hmm_subclusters) {
                    tmp = nchar(reload_info$expected_file_names[[step_count]])
                    tmp = paste0(substr(reload_info$expected_file_names[[step_count]], 1, tmp-13), ".per_chr_subclusters", substr(reload_info$expected_file_names[[step_count]], tmp-12, tmp))
                    saveRDS(subclusters_per_chr, tmp)
                }
            }
            invisible(gc())
            
            if (plot_steps) {
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_tumor_subclusters",step_count),
                         output_filename=sprintf("infercnv.%02d_tumor_subclusters",step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }

            if (inspect_subclusters) {
                plot_subclusters(infercnv_obj,
                                 out_dir=out_dir,
                                 output_filename="infercnv_subclusters")
            }
        }
    }
    else if (analysis_mode != 'subclusters') {
        
        flog.info(sprintf("\n\n\tSTEP %02d: Clustering samples (not defining tumor subclusters)\n", step_count))
        
        if (skip_past < step_count) {
            
            
            infercnv_obj <- define_signif_tumor_subclusters(infercnv_obj,
                                                            p_val=tumor_subcluster_pval,
                                                            hclust_method=hclust_method,
                                                            cluster_by_groups=cluster_by_groups,
                                                            partition_method='none',
                                                            per_chr_hmm_subclusters=FALSE,
                                                            z_score_filter=z_score_filter
                                                            )[[1]]
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
        }
        
    }
    
    
    ## This is a milestone step and results should always be examined here.
    if (skip_past < step_count) {
        if (save_rds) {
            saveRDS(infercnv_obj, file=file.path(out_dir, "preliminary.infercnv_obj"))
        }
        
        invisible(gc())
    
        if (! (no_prelim_plot | no_plot) ) {
            
            #prelim_heatmap_png = "infercnv.preliminary.png"
            
            #if (! file.exists(file.path(out_dir, prelim_heatmap_png))) {
            plot_cnv(infercnv_obj,
                     k_obs_groups=k_obs_groups,
                     cluster_by_groups=cluster_by_groups,
                     cluster_references=cluster_references,
                     plot_chr_scale=plot_chr_scale,
                     chr_lengths=chr_lengths,
                     out_dir=out_dir,
                     title="Preliminary infercnv (pre-noise filtering)",
                     output_filename="infercnv.preliminary", # png ext auto added
                     output_format=output_format,
                     write_expr_matrix=write_expr_matrix,
                     write_phylo=write_phylo,
                     png_res=png_res,
                     useRaster=useRaster)
            #}
        }
    }
    
    ## Below represent optional downstream analysis steps:
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 16
    if (prune_outliers) {
        
        ## ################################
        ## STEP: Remove outliers for viz
        
        flog.info(sprintf("\n\n\tSTEP %02d: Removing outliers\n", step_count))

        if (skip_past < step_count) {
            
            infercnv_obj = remove_outliers_norm(infercnv_obj,
                                                out_method=outlier_method_bound,
                                                lower_bound=outlier_lower_bound,
                                                upper_bound=outlier_upper_bound)
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps) {
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_removed_outliers",step_count),
                         output_filename=sprintf("infercnv.%02d_removed_outliers", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }
        }
    }
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 17
    hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", analysis_mode)
    if (HMM) {

        if (HMM_type == 'i6') {
            hmm_center = 3
            hmm_state_range = c(0,6)
        } else {
            ## i3
            hmm_center = 2
            hmm_state_range = c(1,3)
        }
        
        if (skip_hmm < 1) {
            flog.info(sprintf("\n\n\tSTEP %02d: HMM-based CNV prediction\n", step_count))
            
            if (analysis_mode == 'subclusters') {
                
                if (HMM_type == 'i6') {
                    if (tumor_subcluster_partition_method == "leiden" && per_chr_hmm_subclusters) {
                        hmm.infercnv_obj <- predict_CNV_via_HMM_on_tumor_subclusters_per_chr(infercnv_obj = infercnv_obj, 
                                                                                             subclusters_per_chr = subclusters_per_chr,
                                                                                             t=HMM_transition_prob)
                    }
                    else {
                        hmm.infercnv_obj <- predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj = infercnv_obj,
                                                                                     t=HMM_transition_prob)
                    }
                } else if (HMM_type == 'i3') {
                    hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj,
                                                                                       i3_p_val=HMM_i3_pval,
                                                                                       t=HMM_transition_prob,
                                                                                       use_KS=HMM_i3_use_KS)
                } else {
                    stop("Error, not recognizing HMM_type")
                }

                if (tumor_subcluster_partition_method == 'random_trees') {
                    ## need to redo the hierarchicial clustering, since the subcluster assignments dont always perfectly line up with the top-level dendrogram.
                    hmm.infercnv_obj <- .redo_hierarchical_clustering(hmm.infercnv_obj, hclust_method=hclust_method)
                }
                
            } else if (analysis_mode == 'cells') {
                
                if (HMM_type == 'i6') {
                    hmm.infercnv_obj <- predict_CNV_via_HMM_on_indiv_cells(infercnv_obj, t=HMM_transition_prob)
                } else if (HMM_type == 'i3') {
                    hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_indiv_cells(infercnv_obj,
                                                                                 i3_p_val=HMM_i3_pval,
                                                                                 t=HMM_transition_prob,
                                                                                 use_KS=HMM_i3_use_KS)
                } else {
                    stop("Error, not recognizing HMM_type")
                }
                
                
            } else {
                ## samples mode
                
                if (HMM_type == 'i6') {
                    hmm.infercnv_obj <- predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj,
                                                                                   cluster_by_groups=cluster_by_groups,
                                                                                   t=HMM_transition_prob)
                } else if (HMM_type == 'i3') {
                    hmm.infercnv_obj <- i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj,
                                                                                         cluster_by_groups=cluster_by_groups,
                                                                                         i3_p_val=HMM_i3_pval,
                                                                                         t=HMM_transition_prob,
                                                                                         use_KS=HMM_i3_use_KS
                                                                                         )
                } else {
                    stop("Error, not recognizing HMM_type")
                }

            }
            
            
            ## report predicted cnv regions:
            generate_cnv_region_reports(hmm.infercnv_obj,
                                        output_filename_prefix=sprintf("%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                                        out_dir=out_dir,
                                        ignore_neutral_state=hmm_center,
                                        by=HMM_report_by)

            
            ## ##################################
            ## Note, HMM invercnv object is only leveraged here, but stored as file for future use:
            ## ##################################
            
            if (save_rds) {
                saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
            }

            invisible(gc())
                        
            if (! no_plot) {
                
                ## Plot HMM pred img
                plot_cnv(infercnv_obj=hmm.infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds",step_count),
                         output_filename=sprintf("infercnv.%02d_HMM_pred%s", step_count, hmm_resume_file_token),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=hmm_center,
                         x.range=hmm_state_range,
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
        }
    }

        ## ############################################################
        ## Bayesian Network Mixture Model
        ## ############################################################
        
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 18
    if (skip_hmm < 2) {
        if (HMM == TRUE && BayesMaxPNormal > 0 && length(unique(apply(hmm.infercnv_obj@expr.data,2,unique))) != 1 ) {
            flog.info(sprintf("\n\n\tSTEP %02d: Run Bayesian Network Model on HMM predicted CNVs\n", step_count))
            
            ## the MCMC  object
            
            mcmc_obj <- infercnv::inferCNVBayesNet( infercnv_obj     = infercnv_obj,
                                                   HMM_states        = hmm.infercnv_obj@expr.data,
                                                   file_dir          = out_dir,
                                                   no_plot           = no_plot,
                                                   postMcmcMethod    = "removeCNV",
                                                   out_dir           = file.path(out_dir, sprintf("BayesNetOutput.%s", hmm_resume_file_token)),
                                                   resume_file_token = hmm_resume_file_token,
                                                   quietly           = TRUE,
                                                   CORES             = num_threads,
                                                   plotingProbs      = plot_probabilities,
                                                   diagnostics       = diagnostics,
                                                   HMM_type          = HMM_type,
                                                   k_obs_groups      = k_obs_groups,
                                                   cluster_by_groups = cluster_by_groups,
                                                   reassignCNVs      = reassignCNVs,
                                                   useRaster         = useRaster)

            # mcmc_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj",
            #                                        step_count, hmm_resume_file_token))
            
            if (save_rds) {
                # saveRDS(mcmc_obj, file=mcmc_obj_file)
                saveRDS(mcmc_obj, reload_info$expected_file_names[[step_count]])
            }
        }
    }


        ## ############################################################
        ## Bayesian Filtering
        ## ############################################################

    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 19
    if (skip_hmm < 3) {
        if (HMM == TRUE && BayesMaxPNormal > 0 && length(unique(apply(hmm.infercnv_obj@expr.data,2,unique))) != 1 ) {
            flog.info(sprintf("\n\n\tSTEP %02d: Filter HMM predicted CNVs based on the Bayesian Network Model results and BayesMaxPNormal\n", step_count))
            ## Filter CNV's by posterior Probabilities
            mcmc_obj_hmm_states_list <- infercnv::filterHighPNormals( MCMC_inferCNV_obj = mcmc_obj,
                                                                     HMM_states         = hmm.infercnv_obj@expr.data, 
                                                                     BayesMaxPNormal    = BayesMaxPNormal,
                                                                     useRaster          = useRaster)
            
            hmm_states_highPnormCNVsRemoved.matrix <- mcmc_obj_hmm_states_list[[2]]

            # replace states
            hmm.infercnv_obj@expr.data <- hmm_states_highPnormCNVsRemoved.matrix
            
            ## Save the MCMC inferCNV object
            if (save_rds) {
                saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (! no_plot) {
                ## Plot HMM pred img after cnv removal
                plot_cnv(infercnv_obj=hmm.infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds_Bayes_Net",step_count),
                         output_filename=sprintf("infercnv.%02d_HMM_pred.Bayes_Net.Pnorm_%g",step_count, BayesMaxPNormal),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=hmm_center,
                         x.range=hmm_state_range,
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
            ## write the adjusted CNV report files
            ## report predicted cnv regions:
            adjust_genes_regions_report(mcmc_obj_hmm_states_list[[1]],
                                        # input_filename_prefix=sprintf("%02d_HMM_preds", (step_count-1)),
                                        input_filename_prefix=sprintf("%02d_HMM_pred%s", (step_count-2), hmm_resume_file_token),
                                        output_filename_prefix=sprintf("HMM_CNV_predictions.%s.Pnorm_%g", hmm_resume_file_token, BayesMaxPNormal),
                                        out_dir=out_dir)
        }
    }
        
        
        ## convert from states to representative  intensity values
        
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 20
    if (skip_hmm < 4) {
        if (HMM) {
            flog.info(sprintf("\n\n\tSTEP %02d: Converting HMM-based CNV states to repr expr vals\n", step_count))
            
            if (HMM_type == 'i6') {
                hmm.infercnv_obj <- assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
            } else if (HMM_type == 'i3') {
                hmm.infercnv_obj <- i3HMM_assign_HMM_states_to_proxy_expr_vals(hmm.infercnv_obj)
            }
            
            if (save_rds) {
                saveRDS(hmm.infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            ## Plot HMM pred img
            if (! no_plot) {
                plot_cnv(infercnv_obj=hmm.infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_HMM_preds.repr_intensities",step_count),
                         output_filename=sprintf("infercnv.%02d_HMM_pred%s.Pnorm_%g.repr_intensities", step_count, hmm_resume_file_token, BayesMaxPNormal),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         x.center=1,
                         x.range=c(-1,3),
                         png_res=png_res,
                         useRaster=useRaster
                         )
            }
        }
    }
    
    ## all processes that are alternatives to the HMM prediction wrt DE analysis and/or denoising
    
    ## Step: Filtering significantly DE genes
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 21
    if (mask_nonDE_genes) {
        
        if (!has_reference_cells(infercnv_obj)) {
            stop("Error, cannot mask non-DE genes when there are no normal references set")
        }
        
        flog.info(sprintf("\n\n\tSTEP %02d: Identify and mask non-DE genes\n", step_count))

        if (skip_past < step_count) {

            infercnv_obj <- mask_non_DE_genes_basic(infercnv_obj,
                                                    p_val_thresh=mask_nonDE_pval,
                                                    test.use = test.use,
                                                    center_val=mean(infercnv_obj@expr.data),
                                                    require_DE_all_normals=require_DE_all_normals)
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            
            invisible(gc())
            
            ## Plot incremental steps.
            if (plot_steps) {
                
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         title=sprintf("%02d_mask_nonDE",step_count),
                         output_filename=sprintf("infercnv.%02d_mask_nonDE", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
                
            }
        }
    }
    
    
    if (up_to_step == step_count) {
        flog.info("Reached up_to_step")
        return(infercnv_obj)
    }
    step_count = step_count + 1 # 22
    if (denoise) {
        
        ## ##############################
        ## Step: de-noising
        
        flog.info(sprintf("\n\n\tSTEP %02d: Denoising\n", step_count))

        if (skip_past < step_count) {
            
            if (! is.na(noise_filter)) {
                
                if (noise_filter > 0) {
                    flog.info(paste("::process_data:Remove noise, noise threshold at: ", noise_filter))
                    infercnv_obj <- clear_noise(infercnv_obj,
                                                threshold=noise_filter,
                                                noise_logistic=noise_logistic)
                }
                else {
                    ## noise == 0 or negative...
                    ## don't remove noise.
                }
                
            }
            else {
                ## default, use quantiles, if NA
                flog.info(paste("::process_data:Remove noise, noise threshold defined via ref mean sd_amplifier: ", sd_amplifier))
                infercnv_obj <- clear_noise_via_ref_mean_sd(infercnv_obj,
                                                            sd_amplifier = sd_amplifier,
                                                            noise_logistic=noise_logistic)
            }
            
            if (save_rds) {
                saveRDS(infercnv_obj, reload_info$expected_file_names[[step_count]])
            }
            invisible(gc())
            
            if (plot_steps) {
                plot_cnv(infercnv_obj,
                         k_obs_groups=k_obs_groups,
                         cluster_by_groups=cluster_by_groups,
                         cluster_references=cluster_references,
                         plot_chr_scale=plot_chr_scale,
                         chr_lengths=chr_lengths,
                         out_dir=out_dir,
                         color_safe_pal=FALSE,
                         title=sprintf("%02d_denoised", step_count),
                         output_filename=sprintf("infercnv.%02d_denoised", step_count),
                         output_format=output_format,
                         write_expr_matrix=write_expr_matrix,
                         write_phylo=write_phylo,
                         png_res=png_res,
                         useRaster=useRaster)
            }
            
        }
    }
    
    if (save_final_rds) {
        saveRDS(infercnv_obj, file=file.path(out_dir, "run.final.infercnv_obj"))
    }
    
    if (! no_plot) {
        if (is.null(final_scale_limits)) {
            final_scale_limits = "auto"
        }
        if (is.null(final_center_val)) {
            final_center_val = 1
        }
        
        
        flog.info("\n\n## Making the final infercnv heatmap ##")
        invisible(gc())
        plot_cnv(infercnv_obj,
                 k_obs_groups=k_obs_groups,
                 cluster_by_groups=cluster_by_groups,
                 cluster_references=cluster_references,
                 plot_chr_scale=plot_chr_scale,
                 chr_lengths=chr_lengths,
                 out_dir=out_dir,
                 x.center=final_center_val,
                 x.range=final_scale_limits,
                 title="inferCNV",
                 output_filename="infercnv",
                 output_format=output_format,
                 write_expr_matrix=write_expr_matrix,
                 write_phylo=write_phylo,
                 png_res=png_res,
                 useRaster=useRaster)
    }
    
    return(infercnv_obj)
    
}


#' Subtracting the mean of the reference expr distributions from the observed cells.
#'
#' @title subtract_ref_expr_from_obs()
#'
#' @description Remove the average of the genes of the reference observations from all
#' observations' expression. In the case there are multiple reference groupings,
#' the averages are computed separately for each reference grouping, and the min|max
#' of the averages are subtracted from the observation expression levels. Any values within the range
#' of the min,max of the group are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param inv_log mean values will be determined based on (2^x -1)
#'
#' @param use_bounds if multiple normal data sets are used, it takes the bounds of the means from each set for subtraction.
#'                   Alternatively, will use the mean( mean(normal) for each normal)    default: TRUE
#'
#' @return infercnv_obj containing the reference subtracted values.
#'
#' @keywords internal
#' @noRd
#'

subtract_ref_expr_from_obs <- function(infercnv_obj, inv_log=FALSE, use_bounds=TRUE) {
    ## r = genes, c = cells
    flog.info(sprintf("::subtract_ref_expr_from_obs:Start inv_log=%s, use_bounds=%s", inv_log, use_bounds))
    
    
    if (has_reference_cells(infercnv_obj)) {
        ref_groups = infercnv_obj@reference_grouped_cell_indices
        flog.info("subtracting mean(normal) per gene per cell across all data")
    } else {
        ref_groups = list('proxyNormal' = unlist(infercnv_obj@observation_grouped_cell_indices))
        flog.info("-no reference cells specified... using mean of all cells as proxy")
    }
    
    ref_grp_gene_means <- .get_normal_gene_mean_bounds(infercnv_obj@expr.data, ref_groups, inv_log=inv_log)
    
    infercnv_obj@expr.data <- .subtract_expr(infercnv_obj@expr.data, ref_grp_gene_means, use_bounds=use_bounds)
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- subtract_ref_expr_from_obs(infercnv_obj@.hspike, inv_log=inv_log, use_bounds=use_bounds)
    }
    
    return(infercnv_obj)
    
}

#' @keywords internal
#' @noRd
#'

.get_normal_gene_mean_bounds <- function(expr.data, ref_groups, inv_log=FALSE) {
    
    get_indiv_gene_group_means_bounds_fun <- function(x) {
        
        grp_means = c()
        
        if (inv_log) {
            grp_means <- vapply(ref_groups, function(ref_group) {
                log2(mean(2^x[ref_group] - 1) + 1)
            }, double(1))
        }
        else {
            grp_means <- vapply(ref_groups, function(ref_group) {
                mean(x[ref_group])
            }, double(1))
        }

        names(grp_means) <- names(ref_groups)
        
        return(as.data.frame(t(data.frame(grp_means))))
    }
    
    gene_ref_grp_means <- do.call(rbind, apply(expr.data, 1, get_indiv_gene_group_means_bounds_fun))
    
    rownames(gene_ref_grp_means) <- rownames(expr.data)
    
    return(gene_ref_grp_means)
}


#' @keywords internal
#' @noRd
#'

.subtract_expr <- function(expr_matrix, ref_grp_gene_means, use_bounds=FALSE) {
    
    my.rownames = rownames(expr_matrix)
    my.colnames = colnames(expr_matrix)
    
    flog.info(sprintf("-subtracting expr per gene, use_bounds=%s", use_bounds))
    
    subtract_normal_expr_fun <- function(row_idx) {
        
        gene_means <- as.numeric(ref_grp_gene_means[row_idx, , drop=TRUE])
        
        gene_means_mean <- mean(gene_means)
        
        x <- as.numeric(expr_matrix[row_idx, , drop=TRUE])
        
        row_init = rep(0, length(x))
        
        if (use_bounds) {
            
            grp_min = min(gene_means)
            grp_max = max(gene_means)
            
            above_max = which(x>grp_max)
            below_min = which(x<grp_min)
            
            row_init[above_max] <- x[above_max] - grp_max
            row_init[below_min] <- x[below_min] - grp_min
            ## note, in-between values are left at zero!
            
        } else {
            
            row_init <- x - gene_means_mean
            
        }
        
        return(row_init)
    }
    
    subtr_data <- do.call(rbind, lapply(seq_len(nrow(expr_matrix)), subtract_normal_expr_fun))
    rownames(subtr_data) <- my.rownames
    colnames(subtr_data) <- my.colnames
    
    return(subtr_data)
    
}


#' @description Helper function allowing greater control over the steps in a color palette.
#'              Source: http://menugget.blogspot.com/2011/11/define-color-steps-for-
#'              colorramppalette.html#more
#'
#' @title Helper function allowing greater control over the steps in a color palette.
# 
#' @param steps Vector of colors to change use in the palette
#' @param between Steps where gradients change
#' @param ... Additional arguments of colorRampPalette
#
#' @return Color palette
#'
#' @export
#'
#' @examples
#' color.palette(c("darkblue", "white", "darkred"),
#'               c(2, 2))
#'

color.palette <- function(steps,
                          between=NULL, ...){
    
    if (is.null(between)){
        between <- rep(0, (length(steps) - 1))
    }
    if (length(between) != length(steps) - 1){
        stop("Must have one less \"between\" value than steps")
    }
    
    fill.steps <- cumsum(rep(1, length(steps)) + c(0, between))
    RGB <- matrix(NA, nrow=3, ncol=fill.steps[length(fill.steps)])
    RGB[, fill.steps] <- col2rgb(steps)
    for(i in which(between > 0)){
        col.start <- RGB[, fill.steps[i]]
        col.end <- RGB[, fill.steps[i + 1]]
        for(j in seq(3)){
            vals <- seq(col.start[j],
                        col.end[j],
                        length.out=between[i] + 2)
            vals <- vals[2:(2 + between[i] - 1)]
            RGB[j, (fill.steps[i] + 1):(fill.steps[i + 1] - 1)] <- vals
        }
    }
    new.steps <- rgb(RGB[1, ], RGB[2, ], RGB[3, ], maxColorValue=255)
    pal <- colorRampPalette(new.steps, ...)
    return(pal)
}

# Create a sepList for the heatmap.3 plotting function given integer vectors
# of rows and columns where seperation should take place.
# The expected input to the heatmap function is a list of 2 lists.
# The first list are column based rectangles, and the second row.
# To define a rectangle the index of the row or column where the line of the rectagle
# should be placed is done with a vector of integers, left, bottom, right and top line.
# Ie. list(list(c(1,0,3,10), c(5, 0, 10,10)), list(c(1,2,3,4)))
#
# Args:
# row_count Total number of rows
# col_count Total number of columns
# row_seps Vector of integers indices for row breaks
# col_seps Vector of integer indices for column breaks
#
# Returns
# List of lists of vectors

#' @keywords internal
#' @noRd
#'

create_sep_list <- function(row_count,
                            col_count,
                            row_seps=NULL,
                            col_seps=NULL){
    sepList <- list()
                                        # Heatmap.3 wants a list of boxes for seperating columns
                                        # Column data
    if(!is.null(col_seps) &&
       !any(is.na(col_seps)) &&
       (length(col_seps)>0) &&
       col_count > 0){
        colList <- list()
        for(sep in seq_along(col_seps)){
            colList[[sep]] <- c(col_seps[sep],0,col_seps[sep],row_count)
        }
        sepList[[1]] <- colList
    } else {
        sepList[[1]] <- list()
        sepList[[1]][[1]] <- c(0,0,0,0)
    }
    
                                        # Row data
                                        # This is measured from bottom to top
                                        # So you have to adjust the values of the data
    row_seps <- row_count-row_seps
    if(!is.null(row_seps) &&
       !any(is.na(row_seps)) &&
       (length(row_seps)>0) &&
       row_count > 0){
        rowList <- list()
        for(sep in seq_along(row_seps)){
            rowList[[sep]] <- c(0,row_seps[sep],col_count,row_seps[sep])
        }
        sepList[[2]] <- rowList
    } else {
        sepList[[2]] <- list()
        sepList[[2]][[1]] <- c(0,0,0,0)
    }
    return(sepList)
}


#' @title split_references()
#'
#' @description Split up reference observations in to k groups based on hierarchical clustering.
#'
#'
#' @param infercnv_obj infercnv_object
#'
#' @param num_groups (default: 2)
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

split_references <- function(infercnv_obj,
                             num_groups=2,
                             hclust_method='complete') {
    
    
    flog.info(paste("::split_references:Start", sep=""))
    
    if ("sparseMatrix" %in% is(infercnv_obj@expr.data)) {
        infercnv_obj@expr.data = as.matrix(infercnv_obj@expr.data)
    }

    ref_expr_matrix = infercnv_obj@expr.data[ , get_reference_grouped_cell_indices(infercnv_obj) ]
    
    hc <- hclust(parallelDist(t(ref_expr_matrix), threads=infercnv.env$GLOBAL_NUM_THREADS), method=hclust_method)
    
    split_groups <- cutree(hc, k=num_groups)
    
    ref_groups <- list()
    
    grp_counter = 0
    for (cut_group in unique(split_groups)) {
        grp_counter = grp_counter + 1
        grp_name = sprintf("refgrp-%d", grp_counter)
        cell_names <- names(split_groups[split_groups == cut_group])
        ref_groups[[grp_name]] <- which(colnames(infercnv_obj@expr.data) %in% cell_names)
    }
    
    infercnv_obj@reference_grouped_cell_indices <- ref_groups
    
    return(infercnv_obj)
}


#' @title remove_outliers_norm()
#'
#' @description Set outliers to some upper or lower bound.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param out_method method for computing the outlier bounds (default: "average_bound", involving
#'                   determining the range of values for each cell, and then taking the mean of those bounds.)
#'
#' @param lower_bound setting the lower bound for the data (default: NA, uses out_method above)
#'
#' @param upper_bound setting the upper bound for the data (default: NA, uses out_method above)
#'
#' @return infercnv_obj with data bounds set accordingly.
#'
#' @keywords internal
#' @noRd
#'

remove_outliers_norm <- function(infercnv_obj,
                                 out_method="average_bound",
                                 lower_bound=NA,
                                 upper_bound=NA) {
    
    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))
    
    
    infercnv_obj@expr.data <- .remove_outliers_norm(data=infercnv_obj@expr.data,
                                                    out_method=out_method,
                                                    lower_bound=lower_bound,
                                                    upper_bound=upper_bound)
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- remove_outliers_norm(infercnv_obj@.hspike, out_method, lower_bound, upper_bound)
    }
    
    return(infercnv_obj)
    
}

#' @keywords internal
#' @noRd
#'

.remove_outliers_norm <- function(data,
                                  out_method="average_bound",
                                  lower_bound=NA,
                                  upper_bound=NA) {
    
    flog.info(paste("::remove_outlier_norm:Start",
                    "out_method:", out_method,
                    "lower_bound:" , lower_bound,
                    "upper_bound:", upper_bound))
    
    
    
    if(is.null(data) || nrow(data) < 1 || ncol(data) < 1){
        flog.error("::remove_outlier_norm: Error, something is wrong with the data, either null or no rows or columns")
        stop("Error, something is wrong with the data, either null or no rows or columns")
    }
    
    
    
    if ( (! is.na(lower_bound)) & (! is.na(upper_bound)) ) {
        
        ## using the values as specified.
        flog.info(paste("::remove_outlier_norm: using hard thresholds: ",
                        "lower_bound:" , lower_bound,
                        "upper_bound:", upper_bound) )
        
    } else if (! is.na(out_method)) {
        
                                        # using out_method instead of specified bounds.
        flog.info(paste("::remove_outlier_norm using method:", out_method, "for defining outliers."))
        
        if (out_method == "average_bound"){
            
            bounds = .get_average_bounds(data)
            lower_bound = bounds[1]
            upper_bound = bounds[2]
            
            flog.info(sprintf("outlier bounds defined between: %g - %g", lower_bound, upper_bound))
            
        } else {
            flog.error(paste("::remove_outlier_norm:Error, please",
                             "provide an approved method for outlier",
                             "removal for visualization."))
            stop(991)
        }
    } else {
        flog.error("::remove_outlier_norm:Error, must specify outmethod or define exact bounds")
        stop(992)
    }
    
                                        # apply bounds
    data[data < lower_bound] <- lower_bound
    data[data > upper_bound] <- upper_bound
    
    
    return(data)
}




#' @title center_cell_expr_across_chromosome()
#'
#' @description Centers expression data across all genes for each cell, using the cell mean or median expression
#' value as the center.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param method method to select the center of the cell expression value. (default: 'mean', options: 'mean,'median')
#'
#' @return infercnv_object
#'
#' @keywords internal
#' @noRd
#'

center_cell_expr_across_chromosome <- function(infercnv_obj, method="mean") { # or median
    
    flog.info(paste("::center_smooth across chromosomes per cell"))
    
    infercnv_obj@expr.data <- .center_columns(infercnv_obj@expr.data, method)
    
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- center_cell_expr_across_chromosome(infercnv_obj@.hspike, method)
    }
    
    
    return(infercnv_obj)
}

#' @keywords internal
#' @noRd
#'

.center_columns <- function(expr_data, method) {
    
                                        # Center within columns (cells)
    if (method == "median") {
        row_median <- apply(expr_data, 2, function(x) { median(x, na.rm=TRUE) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_median))
    }
    else {
                                        # by mean
        row_means <- apply(expr_data, 2, function(x) { mean(x, na.rm=TRUE) } )
        
        expr_data <- t(apply(expr_data, 1, "-", row_means))
    }
    return(expr_data)
}




#' @title require_above_min_mean_expr_cutoff ()
#'
#' @description Filters out genes that have fewer than the corresponding mean value across all cell values.
#'
#' @param infercnv_obj  infercnv_object
#'
#' @param min_mean_expr_cutoff  the minimum mean value allowed for a gene to be retained in the expression matrix.
#'
#' @return infercnv_obj  the infercnv_object with lowly or unexpressed genes removed.
#'
#' @keywords internal
#' @noRd
#'

require_above_min_mean_expr_cutoff <- function(infercnv_obj, min_mean_expr_cutoff) {
    
    flog.info(paste("::above_min_mean_expr_cutoff:Start", sep=""))
    
    
    indices <-.below_min_mean_expr_cutoff(infercnv_obj@expr.data, min_mean_expr_cutoff)
    if (length(indices) > 0) {
        flog.info(sprintf("Removing %d genes from matrix as below mean expr threshold: %g",
                          length(indices), min_mean_expr_cutoff))
        
        infercnv_obj <- remove_genes(infercnv_obj, indices)
        
        expr_dim = dim(infercnv_obj@expr.data)
        flog.info(sprintf("There are %d genes and %d cells remaining in the expr matrix.",
                          expr_dim[1], expr_dim[2]))
        
    }
    
    return(infercnv_obj)
    
}

#' @keywords internal
#' @noRd
#'

.below_min_mean_expr_cutoff <- function(expr_data, min_mean_expr) {
    
    average_gene <- rowMeans(expr_data)
    
                                        # Find averages above a certain threshold
    indices <- which(average_gene < min_mean_expr)
    
    return(indices)
    
}




#' @title require_above_min_cells_ref()
#'
#' @description Filters out genes that have fewer than specified number of cells expressing them.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param min_cells_per_gene int indicating number of cells required per gene for both obs and ref data
#'
#' @return infercnv_obj infercnv_object with corresponding genes removed.
#'
#' @keywords internal
#' @noRd
#'

require_above_min_cells_ref <- function(infercnv_obj, min_cells_per_gene) {
    
    genes_passed = which(apply(infercnv_obj@expr.data, 1, function(x) { sum(x>0 & ! is.na(x)) >= min_cells_per_gene}))
    
    num_genes_total = dim(infercnv_obj@expr.data)[1]
    num_removed = num_genes_total - length(genes_passed)
    if (num_removed > 0) {
        
        flog.info(sprintf("Removed %d genes having fewer than %d min cells per gene = %g %% genes removed here",
                          num_removed, min_cells_per_gene, num_removed / num_genes_total * 100))
        
        
        if (num_removed == num_genes_total) {
            
            flog.warn(paste("::All genes removed! Must revisit your data..., cannot continue here."))
            stop(998)
        }
        
        
        infercnv_obj <- remove_genes(infercnv_obj, -1 * genes_passed)
        
        
    }
    else {
        
        flog.info("no genes removed due to min cells/gene filter")
        
    }
    
    return(infercnv_obj)
    
}


#' @title clear_noise()
#'
#' @description Remove values that are too close to the reference cell expr average and are considered noise.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param threshold values within reference mean +- threshold are set to zero.
#'
#' @param noise_logistic uses a logistic (sigmoidal) function to noise removal.
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

clear_noise <- function(infercnv_obj, threshold, noise_logistic=FALSE) {
    
    flog.info(paste("********* ::clear_noise:Start. threshold: ", threshold, sep=""))
    
    if (threshold == 0) {
        return(infercnv_obj); # nothing to do
    }
    
    if (has_reference_cells(infercnv_obj)) {
        ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
        mean_ref_vals = mean(infercnv_obj@expr.data[,ref_idx])
    } else {
        ## no reference
        ## use mean of all data
        mean_ref_vals = mean(infercnv_obj@expr.data)
    }
    
    if (noise_logistic) {
        
        infercnv_obj <- depress_log_signal_midpt_val(infercnv_obj, mean_ref_vals, threshold)
        
    } else {
        
        infercnv_obj@expr.data <- .clear_noise(infercnv_obj@expr.data, threshold, center_pos=mean_ref_vals)
    }
    
    # if (! is.null(infercnv_obj@.hspike)) {
    #     flog.info("-mirroring for hspike")
    #     infercnv_obj@.hspike <- clear_noise(infercnv_obj@.hspike, threshold, noise_logistic)
    # }
    
    return(infercnv_obj)
}

#' @keywords internal
#' @noRd
#'

.clear_noise <- function(expr_data, threshold, center_pos=0) {
    
    upper_bound = center_pos + threshold
    lower_bound = center_pos - threshold
    
    expr_data[expr_data > lower_bound & expr_data < upper_bound] = center_pos
    
    return(expr_data)
}



#' @title clear_noise_via_ref_mean_sd()
#'
#' @description Define noise based on the standard deviation of the reference cell expression data.
#' The range to remove noise would be mean +- sdev * sd_amplifier
#' where sd_amplifier expands the range around the mean to be removed as noise.
#' Data points defined as noise are set to zero.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param sd_amplifier multiplicative factor applied to the standard deviation to alter the noise
#'                     range (default: 1.5)
#'
#' @param noise_logistic uses a logistic (sigmoidal) function to noise removal.
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

clear_noise_via_ref_mean_sd <- function(infercnv_obj, sd_amplifier=1.5, noise_logistic=FALSE) {
    
    if (has_reference_cells(infercnv_obj)) {
        ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
        flog.info("denoising using mean(normal) +- sd_amplifier * sd(normal) per gene per cell across all data")
    }
    else {
        ref_idx = unlist(infercnv_obj@observation_grouped_cell_indices)
        flog.info("-no reference cells specified... using mean and sd of all cells as proxy for denoising")
    }
    vals = infercnv_obj@expr.data[,ref_idx]
    
    mean_ref_vals = mean(vals)
    
    mean_ref_sd <- mean(apply(vals, 2, function(x) sd(x, na.rm=TRUE))) * sd_amplifier
    
    upper_bound = mean_ref_vals + mean_ref_sd
    lower_bound = mean_ref_vals - mean_ref_sd
    
    flog.info(paste(":: **** clear_noise_via_ref_quantiles **** : removing noise between bounds: ",
                    lower_bound, "-", upper_bound, sep=" "))
    
    
    
    if (noise_logistic) {
        
        threshold = mean_ref_sd
        infercnv_obj <- depress_log_signal_midpt_val(infercnv_obj, mean_ref_vals, threshold)
        
    } else {
        # smooth_matrix <- infercnv_obj@expr.data
        
        # smooth_matrix[smooth_matrix > lower_bound & smooth_matrix < upper_bound] = mean_ref_vals
        infercnv_obj@expr.data[infercnv_obj@expr.data > lower_bound & infercnv_obj@expr.data < upper_bound] = mean_ref_vals
        
        # infercnv_obj@expr.data <- smooth_matrix
    }
    
    # if (! is.null(infercnv_obj@.hspike)) {
    #     flog.info("-mirroring for hspike")
    #     infercnv_obj@.hspike <- clear_noise_via_ref_mean_sd(infercnv_obj@.hspike, sd_amplifier, noise_logistic)
    # }
    
    return(infercnv_obj)
}




# Remove the tails of values of a specific chromosome.
# The smooth_matrix values are expected to be in genomic order.
# If the tail is too large and no contig will be left 1/3 of the
# contig is left.
#
# Args:
# smooth_matrix Smoothed values in genomic order.
#                          Row = Genes, Col = Cells.
# chr Indices of the chr in which the tails are to be removed.
# tail_length Length of the tail to remove on both ends of the
#                        chr indices.
# Returns:
# Indices to remove.

#' @keywords internal
#' @noRd
#'


.remove_tails <- function(smooth_matrix, chr, tail_length){
    
                                        #flog.info(paste("::remove_tails:Start.", sep=""))
    chr_length <- length(chr)
    if ((tail_length < 3) || (chr_length < 3)){
        return(c())
    }
    if (chr_length < (tail_length * 2)){
        tail_length <- floor(chr_length / 3)
    }
    remove_indices <- chr[seq_len(tail_length)]
    remove_indices <- c(remove_indices,
                        chr[ ( (chr_length + 1) - tail_length):
                             chr_length])
    
    return(remove_indices)
}


#' @title smooth_by_chromosome()
#'
#' @description Smooth expression values for each cell across each chromosome by using a
#' moving average with a window of specified length.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of window (number of genes) for the moving average
#'
#' @param smooth_ends perform smoothing at the ends of the chromosomes (default: TRUE)
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

smooth_by_chromosome <- function(infercnv_obj, window_length, smooth_ends=TRUE) {
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window(data=chr_data,
                                               window_length=window_length)
            
            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))
            
            infercnv_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- smooth_by_chromosome(infercnv_obj@.hspike, window_length, smooth_ends)
    }
    
    
    return(infercnv_obj)
}

#' @keywords internal
#' @noRd
#'

.smooth_window <- function(data, window_length) {
    
    flog.debug(paste("::smooth_window:Start.", sep=""))
    
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }
    
    ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
    
    data_sm <- apply(data,
                     2,
                     .smooth_helper,
                     window_length=window_length)
    
    
    ## Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    
    
    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))
    
    
    return(data_sm)
}


# Helper function for smoothing with a moving average.
#
# Args:
# obs_data: Data to smooth
# window_length:  Length of the window to smooth.
#
# Returns:
# Data smoothed.
##.smooth_helper <- function(obs_data, tail_length) {

#' @keywords internal
#' @noRd
#'

.smooth_helper <- function(obs_data, window_length) {
                                        # strip NAs out and replace after smoothing
    orig_obs_data = obs_data
    
    nas = is.na(obs_data)
    
    obs_data = obs_data[!nas]
    
    obs_length <- length(obs_data)
    end_data <- obs_data
    
    tail_length = (window_length - 1)/2
    if (obs_length >= window_length) {
        end_data <- .smooth_center_helper(obs_data, window_length)
    }
    
                                        # end_data will have the end positions replaced with mean values, smoothing just at the ends.
    
    obs_count <- length(obs_data)
    
    numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))
    
                                        # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
    iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))
    
    for (tail_end in seq_len(iteration_range)) {
        end_tail = obs_count - tail_end + 1
        
        d_left = tail_end - 1
        d_right = obs_count - tail_end
        d_right = ifelse(d_right > tail_length, tail_length, d_right)
        
        r_left = tail_length - d_left
        r_right = tail_length - d_right
        
        denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)
        
        left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
        right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]
        
        numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
        
        end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
        end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
    }
    
    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data
    
    return(orig_obs_data)
}

smooth_by_chromosome_coordinates <- function(infercnv_obj, window_length) {
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            smoothed_chr_data = .smooth_window_by_coordinates(data=chr_data,
                                                              window_length=window_length,
                                                              genes=infercnv_obj@gene_order[which(infercnv_obj@gene_order$chr == chr) ,])
            
            flog.debug(paste0("smoothed data: ", paste(dim(smoothed_chr_data), collapse=",")))
            
            infercnv_obj@expr.data[chr_genes_indices, ] <- smoothed_chr_data
        }
    }
    
    if (! is.null(infercnv_obj@.hspike)) {  ## for now not supported.
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- smooth_by_chromosome_coordinates(infercnv_obj@.hspike, window_length=51)
    }
    
    
    return(infercnv_obj)
}

.smooth_window_by_coordinates <- function(data, window_length, genes) {
    
    flog.debug(paste("::smooth_window:Start.", sep=""))
    
    if (window_length < 2){
        flog.warn("window length < 2, returning original unmodified data")
        return(data)
    }
    
    ## Fix ends that couldn't be smoothed since not spanned by win/2 at ends.
    
    data_sm <- apply(data,
                     2,
                     .smooth_helper_by_coordinates,
                     window_length=window_length,
                     genes=genes)
    
    
    ## Set back row and column names
    row.names(data_sm) <- row.names(data)
    colnames(data_sm) <- colnames(data)
    
    
    flog.debug(paste("::smooth_window: dim data_sm: ", dim(data_sm), sep=" "))
    
    
    return(data_sm)
}

.smooth_helper_by_coordinates <- function(obs_data, window_length=10000000, genes) {

    ### No handling of NAs for this method of smoothing.

    end_data <- obs_data

    for (i in seq_along(obs_data)) {

        current_pos = (genes$start[i] + genes$stop[i]) / 2

        around_indices = which((genes$start > (current_pos - window_length)) & (genes$stop < (current_pos + window_length)))
        if (length(around_indices) == 0) {  ## in case the window is smaller than the current gene's size
            around_indices = i
        }
        weights = 1 - abs(((genes$stop[around_indices] + genes$start[around_indices])/2) - current_pos)/window_length
        if (length(around_indices < 10)) {
            to_add = floor(length(around_indices - 10)/2)

            new_low = max(1, (min(around_indices) - to_add))
            new_high = min(length(obs_data), (max(around_indices) + to_add))
            weights = c(rep(0.1, (min(around_indices) - new_low)), weights, rep(0.1, (new_high - max(around_indices))))
            around_indices = seq(new_low, new_high)
        }

        end_data[i] = sum(obs_data[around_indices] * weights) / sum(weights)
    }

    return(end_data)
}



# Smooth vector of values over the given window length.
#
# Args:
# obs_data Vector of data to smooth with a moving average.
# window_length Length of the window for smoothing.
#        Must be and odd, positive, integer.
#
# Returns:
# Vector of values smoothed with a moving average.

#' @keywords internal
#' @noRd
#'

.smooth_center_helper <- function(obs_data, window_length){
    
    nas = is.na(obs_data)
    vals = obs_data[! nas]
    
    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))
    
    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)
    
                                        #    flog.info(paste("custom filter = ", custom_filter, "\n and window_length =", window_length, "\n and nrow(data) = ", nrow(obs_data), sep=""))
    
                                        #smoothed = filter(vals, rep(1 / window_length, window_length), sides=2)
    smoothed = filter(vals, custom_filter, sides=2)
    
    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]
    
    obs_data[! nas] = vals
    
    return(obs_data)
}


#' @title smooth_by_chromosome_runmeans
#'
#' @description uses the simpler caTools:runmeans() to perform smoothing operations.
#'
#' @param infercnv_obj infercnv object
#'
#' @param window_length window length to use for smoothing.
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'


smooth_by_chromosome_runmeans <- function(infercnv_obj, window_length) {
    
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    
    for (chr in chrs) {
        chr_genes_indices = which(gene_chr_listing == chr)
        flog.info(paste0("smooth_by_chromosome: chr: ",chr))
        
        chr_data=infercnv_obj@expr.data[chr_genes_indices, , drop=FALSE]
        
        if (nrow(chr_data) > 1) {
            chr_data = apply(chr_data, 2, caTools::runmean, k=window_length)
            
            infercnv_obj@expr.data[chr_genes_indices, ] <- chr_data
        }
    }
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- smooth_by_chromosome_runmeans(infercnv_obj@.hspike, window_length)
    }
    
    
    return(infercnv_obj)
}






#' @title get_average_bounds()
#'
#' @description Computes the mean of the upper and lower bound for the data across all cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return (lower_bound, upper_bound)
#'
#' @keywords internal
#' @noRd
#'

get_average_bounds <- function (infercnv_obj) {
    
    return(.get_average_bounds(infercnv_obj@expr.data))
    
}


#' @keywords internal
#' @noRd
#'

.get_average_bounds <- function(expr_matrix) {
    
    lower_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[1]]))
    upper_bound <- mean(apply(expr_matrix, 2,
                              function(x) quantile(x, na.rm=TRUE)[[5]]))
    
    return(c(lower_bound, upper_bound))
}

#' @title log2xplus1()
#'
#' @description Computes log(x+1), updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

log2xplus1 <- function(infercnv_obj) {
    
    flog.info("transforming log2xplus1()")
    
    infercnv_obj@expr.data <- log2(infercnv_obj@expr.data + 1)
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- log2xplus1(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
    
}



#' @title invert_log2xplus1()
#'
#' @description Computes 2^x - 1
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

invert_log2xplus1 <- function(infercnv_obj) {
    
    flog.info("inverting log2xplus1()")
    
    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data - 1
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- invert_log2xplus1(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
}


#' @title invert_log2()
#'
#' @description Computes 2^x
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

invert_log2 <- function(infercnv_obj) {
    
    flog.info("invert_log2(), computing 2^x")
    
    infercnv_obj@expr.data <- 2^infercnv_obj@expr.data
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- invert_log2(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
}


#' @title make_zero_NA()
#'
#' @description Converts zero to NA
#' Updates infercnv_obj@expr.data
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

make_zero_NA <- function(infercnv_obj) {
    
    flog.info("make_zero_NA()")
    
    infercnv_obj@expr.data <- infercnv_obj@expr.data[infercnv_obj@expr.data == 0] <- NA
    
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- make_zero_NA(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
    
}


#' @title transform_to_reference_based_Zscores()
#'
#' @description Computes mean and standard deviation for the reference cells, then uses these values
#' to compute gene Z-scores for both the reference and observation cells.
#'
#' Note, reference cell gene expression values will be centered at zero after this operation.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

transform_to_reference_based_Zscores <- function(infercnv_obj) {

    ## center and convert to z-scores
    flog.info(paste("::center_and_Zscore_conversion", sep=""))
    
    ## remember, genes are rows, cells are cols
    
    ## centering and z-scores based on the reference (normal) cells:

    ## ref data represent the null distribution
    ref_idx = get_reference_grouped_cell_indices(infercnv_obj)
    
    ref_data = infercnv_obj@expr.data[,ref_idx]
    
    gene_ref_mean = apply(ref_data, 1, function(x) {mean(x, na.rm=TRUE)})
    gene_ref_sd = apply(ref_data, 1, function(x) {sd(x, na.rm=TRUE)})
    
    ## assume at least Poisson level variation
    gene_ref_sd = pmax(gene_ref_sd, gene_ref_mean)
    
    ## center all genes at the ref (normal) center:
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_mean, FUN="-")
    
    ## convert to z-scores based on the ref (normal) distribution
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 1, gene_ref_sd, FUN="/") # make all data z-scores based on the ref data distribution.
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- transform_to_reference_based_Zscores(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
    
}


#' @title mean_center_gene_expr()
#'
#' @description mean-center all gene expression values across all cells
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

mean_center_gene_expr <- function(infercnv_obj) {
    
    flog.info(paste("::centering", sep=""))
    
    infercnv_obj@expr.data <- sweep(infercnv_obj@expr.data, 1, rowMeans(infercnv_obj@expr.data, na.rm=TRUE), FUN="-")
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- mean_center_gene_expr(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
}


#' @title get_reference_grouped_cell_indices()
#'
#' @description Retrieves the matrix indices for the columns correspoinding to the reference cells.
#'
#' @param infercnv_obj infercnv_object
#'
#' @return vector of column indices
#'
#' @keywords internal
#' @noRd
#'

get_reference_grouped_cell_indices <- function(infercnv_obj) {
    
    return( unlist(infercnv_obj@reference_grouped_cell_indices) )
    
}


#' @title apply_max_threshold_bounds()
#'
#' @description Assumes centered at zero and sets bounds to +- threshold value.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param threshold value to threshold the data
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

apply_max_threshold_bounds <- function(infercnv_obj, threshold) {
    
    flog.info(paste("::process_data:setting max centered expr, threshold set to: +/-: ", threshold))
    
    infercnv_obj@expr.data[infercnv_obj@expr.data > threshold] <- threshold
    infercnv_obj@expr.data[infercnv_obj@expr.data < (-1 * threshold)] <- -1 * threshold
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- apply_max_threshold_bounds(infercnv_obj@.hspike, threshold)
    }
    
    return(infercnv_obj)
}


#' @title remove_genes_at_ends_of_chromosomes()
#'
#' @description Removes genes that are within window_length/2 of the ends of each chromosome.
#'
#' @param infercnv_obj infercnv_object
#'
#' @param window_length length of the window to use.
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

remove_genes_at_ends_of_chromosomes <- function(infercnv_obj, window_length) {
    
    contig_tail = (window_length - 1) / 2
    
    remove_indices <- c()
    gene_chr_listing = infercnv_obj@gene_order[[C_CHR]]
    chrs = unlist(unique(gene_chr_listing))
    for (chr in chrs){
                                        #flog.info(paste("::process_data:Remove tail contig ",chr, ".", sep=""))
        remove_chr <- .remove_tails(infercnv_obj@expr.data,
                                    which(gene_chr_listing == chr),
                                    contig_tail)
        
                                        #flog.debug(paste("::process_data:Remove tail - removing indices for chr: ", chr, ", count: ", length(remove_chr), sep=""))
        
        remove_indices <- c(remove_indices, remove_chr)
        
    }
    if (length(remove_indices) > 0){
        
        infercnv_obj = remove_genes(infercnv_obj, remove_indices)
        
        flog.info(paste("::process_data:Remove genes at chr ends, ",
                        "new dimensions (r,c) = ",
                        paste(dim(infercnv_obj@expr.data), collapse=","),
                        " Total=", sum(infercnv_obj@expr.data, na.rm=TRUE),
                        " Min=", min(infercnv_obj@expr.data, na.rm=TRUE),
                        " Max=", max(infercnv_obj@expr.data, na.rm=TRUE),
                        ".", sep=""))
        
    }
    else {
        flog.error("No genes removed at chr ends.... something wrong here")
        stop(1234)
    }
    
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- remove_genes_at_ends_of_chromosomes(infercnv_obj@.hspike, window_length)
    }
    
    return(infercnv_obj)
    
}


#' @title normalize_counts_by_seq_depth()
#'
#' @description Normalizes count data by total sum scaling
#'
#' For single cell data, a typical normalization factor is 1e5, providing counts per 100k total counts.
#' If a normalization factor is not provided, the median lib size is used.:
#'
#' @param infercnv_obj infercnv_object
#'
#' @param normalize_factor  total counts to scale the normalization to (default: NA, computed as described above)
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

normalize_counts_by_seq_depth <- function(infercnv_obj, normalize_factor=NA) {
    
    data <- infercnv_obj@expr.data
    
    normalized_data <- .normalize_data_matrix_by_seq_depth(data, normalize_factor)
    
    
    infercnv_obj@expr.data <- normalized_data
    
    return(infercnv_obj)
    
}


#' @keywords internal
#' @noRd
#'

.normalize_data_matrix_by_seq_depth <- function(counts.matrix, normalize_factor=NA) {
    
    flog.info("normalizing counts matrix by depth")
    
    data <- counts.matrix
    
    cs = colSums(data)
    
    ## make fraction of total counts:
    data <- sweep(data, STATS=cs, MARGIN=2, FUN="/")
    
    if (is.na(normalize_factor)) {

        normalize_factor = median(cs)
        
        flog.info(sprintf("Computed total sum normalization factor as median libsize: %f", normalize_factor))
        
    } else {
        flog.info(sprintf("Using specified normalization factor: %f", normalize_factor))
    }
    
    if (is.na(normalize_factor)) {
        stop("Error, normalize factor not estimated")
    }
    
    data <- data * normalize_factor
    
    return(data)
    
}



#' @title anscombe_transform()
#'
#' @description Performs Anscombe's transformation:
#'    y = 2 * sqrt(x + 3/8)
#' as per
#' https://en.wikipedia.org/wiki/Anscombe_transform
#'
#' @param infercnv_obj infercnv_object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

anscombe_transform <- function(infercnv_obj) {
    
    infercnv_obj@expr.data <- 2 * sqrt(infercnv_obj@expr.data + 3/8)
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- anscombe_transform(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
    
}

#' @keywords internal
#' @noRd
#'
add_pseudocount <- function(infercnv_obj, pseudocount) {
    
    flog.info(sprintf("Adding pseudocount: %g", pseudocount))
    
    infercnv_obj@expr.data = infercnv_obj@expr.data + pseudocount
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- add_pseudocount(infercnv_obj@.hspike, pseudocount)
    }
    
    return(infercnv_obj)
}


#' @title scale_infercnv_expr
#'
#' @description performs scaling to expression values for each cell,
#' assigning all values to a standard normal centered at zero.
#'
#' @param infercnv_obj infercnv object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

scale_infercnv_expr <- function(infercnv_obj) {
    
    flog.info("-scaling expr data")
    infercnv_obj@expr.data = t(scale(t(infercnv_obj@expr.data)))
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- scale_infercnv_expr(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
}



#' @keywords internal
#' @noRd
#'

cross_cell_normalize <- function(infercnv_obj) {
    
    ## using upper quartile normalization
    
    flog.info("-cross cell normalization")
    
    upper_quart = apply(infercnv_obj@expr.data, 2, quantile, probs=0.75)
    mean_upper_quart = mean(upper_quart)
    infercnv_obj@expr.data = sweep(infercnv_obj@expr.data, 2, mean_upper_quart/upper_quart, "*")
    
    
    if (! is.null(infercnv_obj@.hspike)) {
        flog.info("-mirroring for hspike")
        infercnv_obj@.hspike <- cross_cell_normalize(infercnv_obj@.hspike)
    }
    
    return(infercnv_obj)
    
    
}



#' @title .redo_hierarchical_clustering()
#'
#' @description Recomputes hierarchical clustering for subclusters
#'
#'
#' @param infercnv_obj infercnv_object
#'
#' @param hclust_method clustering method to use (default: 'complete')
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'

.redo_hierarchical_clustering <- function(infercnv_obj, hclust_method) {

    subclusters = infercnv_obj@tumor_subclusters$subclusters

    for (grp_name in names(subclusters)) {

        grp_cell_idx = unlist(subclusters[[grp_name]], recursive=TRUE)
        grp_cell_idx = grp_cell_idx[order(grp_cell_idx)] # retain relative ordering in the expr matrix
        
        grp_expr_data = infercnv_obj@expr.data[, grp_cell_idx, drop=FALSE]
        
        hc <- hclust(parallelDist(t(grp_expr_data), threads=infercnv.env$GLOBAL_NUM_THREADS), method=hclust_method)

        infercnv_obj@tumor_subclusters$hc[[grp_name]] <- hc
    }

    return(infercnv_obj)
    
}


# This function returns TRUE wherever elements are the same, including NA's,
# and FALSE everywhere else.
compareNA <- function(v1,v2) {
    if ((is.null(v1) != is.null(v2)) || (length(v1) != length(v2))) {
        return(FALSE)
    }
    if (all((v1 == v2) | (is.na(v1) & is.na(v2)))) {
        return(TRUE)
    }
    return(FALSE)
}

#' @return FALSE if relevant args are not identical, TRUE if they are
#'
#' @keywords internal
#' @noRd
#'

.compare_args <- function(current_args, relevant_args, loaded_options) {

    diff_a = any(relevant_args %in% names(current_args[!(names(current_args) %in% names(loaded_options))]))
    diff_b = any(relevant_args %in% names(loaded_options[!(names(loaded_options) %in% names(current_args))]))

    if (diff_a || diff_b) {
        return(FALSE)
    }

    return(all(sapply(relevant_args[relevant_args %in% names(current_args)], function(arg_n) {
        all(compareNA(current_args[[arg_n]], loaded_options[[arg_n]]))
    })))
}


#' @keywords internal
#' @noRd
#'

.get_relevant_args_list <- function(
# .get_relevant_args_list <- function(out_dir=NULL,
#                                     HMM=FALSE,
#                                     HMM_type='i6',
#                                     num_ref_groups=NULL,
#                                     analysis_mode='samples',
#                                     tumor_subcluster_partition_method,
#                                     smooth_method='pyramidinal',
#                                     max_centered_threshold=3,
#                                     remove_genes_at_chr_ends=FALSE,
#                                     prune_outliers=FALSE,
#                                     BayesMaxPNormal=0.5,
#                                     mask_nonDE_genes=FALSE,
#                                     denoise=FALSE,
#                                     noise_filter=NA,
#                                     sd_amplifier=1.5,
#                                     noise_logistic=FALSE,
                                    ...) {
    run_arguments = list(...)

    # creation args  ## add check for matrix size, col/row names, and ref/obs indices
    relevant_args = vector("list", 21) # 21 = max steps
    expected_file_names = vector("character", 21)
    step_i = 1

    # 1 _incoming_data.infercnv_obj
    relevant_args[[step_i]] = c("chr_exclude", "max_cells_per_group", "min_max_counts_per_cell", "counts_md5")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_incoming_data.infercnv_obj", step_i))
    step_i = step_i + 1

    # 2 _reduced_by_cutoff.infercnv_obj
    relevant_args[[step_i]] = c("cutoff", "min_cells_per_gene")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_reduced_by_cutoff.infercnv_obj", step_i))
    step_i = step_i + 1

    # 3 _normalized_by_depth%s.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type") 
        if (run_arguments$HMM_type == "i6") {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "sim_method", "hspike_aggregate_normals", "sim_foreground")
        }
    }
    resume_file_token = ifelse((run_arguments$HMM), paste0("HMM", run_arguments$HMM_type), "")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_normalized_by_depth%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 4 _logtransformed%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_logtransformed%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 5 _scaled%s.infercnv_obj
    relevant_args[[step_i]] = c("scale_data")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_scaled%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 6 _split_%sf_refs%s.infercnv_obj
    relevant_args[[step_i]] = c("num_ref_groups")
    if (!is.null(run_arguments$num_ref_groups)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_split_%sf_refs%s.infercnv_obj", step_i, resume_file_token, run_arguments$num_ref_groups))
    }
    step_i = step_i + 1

    # 7 _tumor_subclusters%s.%s.infercnv_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method == 'random_trees') {
        resume_file_token = paste0(resume_file_token, ".rand_trees")
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method", "tumor_subcluster_pval", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.%s.infercnv_obj", step_i, resume_file_token, run_arguments$tumor_subcluster_partition_method))
    }
    step_i = step_i + 1

    # 8 _remove_ref_avg_from_obs_logFC%s.infercnv_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_logFC%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 9 _apply_max_centered_expr_threshold%s.infercnv_obj
    relevant_args[[step_i]] = c("max_centered_threshold")
    if(!is.na(run_arguments$max_centered_threshold)) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "max_centered_threshold")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_apply_max_centered_expr_threshold%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 10 _smoothed_by_chr%s.infercnv_obj
    relevant_args[[step_i]] = c("smooth_method", "window_length")
    if (run_arguments$smooth_method == 'pyramidinal') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "smooth_ends")
    }
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_smoothed_by_chr%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 11 _recentered_cells_by_chr%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_recentered_cells_by_chr%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 12 _remove_ref_avg_from_obs_adjust%s.infercnv_obj
    relevant_args[[step_i]] = c("ref_subtract_use_mean_bounds")
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_ref_avg_from_obs_adjust%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 13 _remove_gene_at_chr_ends%s.infercnv_obj
    relevant_args[[step_i]] = c("remove_genes_at_chr_ends")
    if (run_arguments$remove_genes_at_chr_ends == TRUE) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "window_length")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_gene_at_chr_ends%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 14 _invert_log_transform%s.infercnv_obj
    relevant_args[[step_i]] = c()
    expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_invert_log_transform%s.infercnv_obj", step_i, resume_file_token))
    step_i = step_i + 1

    # 15 _tumor_subclusters%s.infercnv_obj
    relevant_args[[step_i]] = c("analysis_mode", "tumor_subcluster_partition_method", "z_score_filter")
    if (run_arguments$analysis_mode == 'subclusters' & run_arguments$tumor_subcluster_partition_method != 'random_trees') {
        resume_file_token = paste0(resume_file_token, '.', run_arguments$tumor_subcluster_partition_method)
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "k_nn", "leiden_resolution", "tumor_subcluster_pval", "leiden_method", "leiden_function", "leiden_method_per_chr", "leiden_function_per_chr", "leiden_resolution_per_chr", "per_chr_hmm_subclusters", "per_chr_hmm_subclusters_references", "hclust_method", "cluster_by_groups") 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_tumor_subclusters%s.infercnv_obj", step_i, resume_file_token))
    }
    else if (run_arguments$analysis_mode != 'subclusters') {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "tumor_subcluster_pval", "hclust_method", "cluster_by_groups")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_no_subclustering%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 16 _remove_outlier%s.infercnv_obj
    relevant_args[[step_i]] = c("prune_outliers")
    if (run_arguments$prune_outliers) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "outlier_method_bound", "outlier_lower_bound", "outlier_upper_bound")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_remove_outlier%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 17 _HMM_pred%s.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_type", "analysis_mode", "HMM_transition_prob", "HMM_report_by")
        hmm_resume_file_token = paste0(resume_file_token, ".hmm_mode-", run_arguments$analysis_mode)
        # hmm.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj", step_i, hmm_resume_file_token))  ########
        if (run_arguments$analysis_mode == 'subclusters' && run_arguments$tumor_subcluster_partition_method == 'random_trees') {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "hclust_method")
        }

        if (run_arguments$HMM_type == "i6") {
            # relevant_args[[step_i]] = c(relevant_args[[step_i]], )

        }
        else {
            relevant_args[[step_i]] = c(relevant_args[[step_i]], "HMM_i3_pval", "HMM_i3_use_KS")
        }
        expected_file_names[[step_i]] = hmm.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred%s.infercnv_obj", step_i, hmm_resume_file_token))
    }
    step_i = step_i + 1

    # 18 _HMM_pred.Bayes_Net%s.mcmc_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
    }
    step_i = step_i + 1

    # 19 _HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj
    relevant_args[[step_i]] = c("HMM", "BayesMaxPNormal")
    if (run_arguments$HMM == TRUE & run_arguments$BayesMaxPNormal > 0) {
        # relevant_args[[step_i]] = c(relevant_args[[step_i]], "diagnostics", "reassignCNVs")
        # mcmc_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token))
        # mcmc.infercnv_obj_file = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        # expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.mcmc_obj", step_i, hmm_resume_file_token)), 
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.Bayes_Net%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 20 _HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj
    relevant_args[[step_i]] = c("HMM")
    if (run_arguments$HMM) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "BayesMaxPNormal")
        # hmm.infercnv_obj_file = file.path(out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, BayesMaxPNormal))
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_HMM_pred.repr_intensities%s.Pnorm_%g.infercnv_obj", step_i, hmm_resume_file_token, run_arguments$BayesMaxPNormal))
    }
    step_i = step_i + 1

    # 21 _mask_nonDE%s.infercnv_obj
    relevant_args[[step_i]] = c("mask_nonDE_genes")
    if (run_arguments$mask_nonDE_genes) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "require_DE_all_normals", "test.use", "mask_nonDE_pval")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_mask_nonDE%s.infercnv_obj", step_i, resume_file_token))
    }
    step_i = step_i + 1

    # 22 _denoise%s.NF_%s.SD_%g.NL_%s.infercnv_obj
    relevant_args[[step_i]] = c("denoise")
    if (run_arguments$denoise) {
        relevant_args[[step_i]] = c(relevant_args[[step_i]], "noise_filter", "sd_amplifier", "noise_logistic")
        expected_file_names[[step_i]] = file.path(run_arguments$out_dir, sprintf("%02d_denoise%s.NF_%s.SD_%g.NL_%s.infercnv_obj", step_i, resume_file_token, run_arguments$noise_filter, run_arguments$sd_amplifier, run_arguments$noise_logistic))
    }
    step_i = step_i + 1

    return(list(relevant_args=relevant_args, expected_file_names=expected_file_names))
}

