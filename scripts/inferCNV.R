#!/usr/bin/env Rscript

# To use inferCNV via command-line interface, first install inferCNV per Wiki,
# then run a command like the following:
#
# ./inferCNV.R \
#     --raw_counts_matrix="../example/oligodendroglioma_expression_downsampled.counts.matrix" \
#     --annotations_file="../example/oligodendroglioma_annotations_downsampled.txt" \
#     --gene_order_file="../example/gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt" \
#     --ref_group_names="Microglia/Macrophage,Oligodendrocytes (non-malignant)" \
#     --cutoff=1 \
#     --out_dir="output_cli" \
#     --cluster_by_groups \
#     --denoise
#     --median_filter

# Load libraries
library(optparse)
library(futile.logger)
#if (!require('fastcluster')) {
#    warning("fastcluster library not available, using the default hclust method instead.")
#}
library(infercnv)
options("preferRaster" = TRUE)
options(scipen = 100)

# Logging level choices
C_LEVEL_CHOICES <- c("DEBUG", "INFO", "WARN", "ERROR")
# Visualization outlier thresholding and bounding method choices
C_VIS_OUTLIER_CHOICES <- c("average_bound")
C_REF_SUBTRACT_METHODS <- c("by_mean", "by_quantiles")
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")

CHR = "chr"
START = "start"
STOP = "stop"


flog.threshold("INFO") #initialize to info setting.

#' Check arguments and make sure the user input meet certain
#' additional requirements.
#'
#' Args:
#'    @param arguments: Parsed arguments from user.
#'
#' Returns:
#'    @return: Updated arguments
check_arguments <- function(arguments){

    flog.info(paste("::check_arguments:Start", sep=""))
    # Require the name of a output pdf file
    if ( (!( "output_dir" %in% names(arguments))) || (arguments$output_dir == "") || (is.na(arguments$output_dir)) ) {
        flog.error(paste(":: --output_dir: Please enter a file path to ",
                                "save the heatmap.",
                                 sep=""))

        stop("error, no --output_dir")
    }

    # Require the cut off to be above 0
    if (arguments$cutoff < 0){
        flog.error(paste(":: --cutoff: Please enter a value",
                                "greater or equal to zero for the cut off.",
                                sep=""))

        stop("error, no --cutoff")
    }

    # Require the logging level to be one handled by logging
    if (!(arguments$log_level %in% C_LEVEL_CHOICES)){
        flog.error(paste(":: --log_level: Please use a logging level ",
                                "given here: ", C_LEVEL_CHOICES,
                                collapse=",", sep=""))
        stop("error, not recognizing log level")
    }
    flog.threshold(arguments$log_level)

    # Require the visualization outlier detection to be a correct choice.
    if (!(arguments$bound_method_vis %in% C_VIS_OUTLIER_CHOICES)){
        flog.error(paste(":: --vis_bound_method: Please use a method ",
                                "given here: ", C_VIS_OUTLIER_CHOICES,
                                collapse=",", sep=""))
        stop("error, must specify acceptable --vis_bound_method")
    }

    if (! (arguments$ref_subtract_method %in% C_REF_SUBTRACT_METHODS) ) {
        flog.error(paste(":: --ref_subtract_method: acceptable values are: ",
                                paste(C_REF_SUBTRACT_METHODS, collapse=","), sep="") )
        stop("error, must specify acceptable --ref_subtract_method")
    }


    if (! (arguments$hclust_method %in% C_HCLUST_METHODS) ) {
        flog.error(paste(":: --hclust_method: acceptable values are: ",
                                paste(C_HCLUST_METHODS, collapse=","), sep="") )
        stop("error, must specify acceptable --hclust_method")
    }

    # # Warn that an average of the samples is used in the absence of
    # # normal / reference samples
    # if (is.null(arguments$reference_observations)){
    #     logging::logwarn(paste(":: --reference_observations: No reference ",
    #                   "samples were given, the average of the samples ",
    #                   "will be used.",
    #                   sep=""))
    # }

    # Make sure the threshold is centered.
    arguments$max_centered_expression <- abs(arguments$max_centered_expression)
    arguments$magnitude_filter <- abs(arguments$magnitude_filter)

    # Require the contig tail to be above 0
    if (is.na(arguments$contig_tail)){
        arguments$contig_tail <- (arguments$window_length - 1) / 2
    }

    if (arguments$contig_tail < 0){
        flog.error(paste(":: --tail: Please enter a value",
                                "greater or equal to zero for the tail.",
                                sep=""))

        stop(980)
    }

    # if (! is.na(suppressWarnings(as.integer(arguments$name_ref_groups)))){
    #     arguments$name_ref_groups <- list(as.integer(arguments$name_ref_groups))
    # } else {

    if (! is.na(arguments$name_ref_groups)) {
        arguments$name_ref_groups <- unlist(strsplit(arguments$name_ref_groups,","))
    }
    else {
        if(!is.null(arguments$num_ref_groups)) {
            flog.error(paste("::  cannot use --num_ref_groups without",
                                    "using --ref_groups as the average of all",
                                    "cells is used."))
            stop(978)
        }
    }
    return(arguments)
}

# Command line arguments
pargs <- optparse::OptionParser(usage=paste("%prog [options]",
                                            "--output_dir directory",
                                            "data_matrix genomic_positions"))

pargs <- optparse::add_option(pargs, c("--cutoff"),
                              type="numeric",
                              default=1,
                              action="store",
                              dest="cutoff",
                              metavar="Cutoff",
                              help=paste("A number >= 0 is expected. Cut-off for the min",
                                         "average read counts per gene among reference cells.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--min_cells_per_gene"),
                              type="numeric",
                              default=3,
                              action="store",
                              dest="min_cells_per_gene",
                              metavar="Minimum cells per gene",
                              help=paste("minimum number of reference cells requiring ",
                                         "expression measurements to include the ",
                                         "corresponding gene. [Default %default]"))

pargs <- optparse::add_option(pargs, c("--out_dir"),
                              type="character",
                              default=".",
                              action="store",
                              dest="out_dir",
                              metavar="Output Directory",
                              help=paste("Path to directory to deposit outputs. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--window_length"),
                              type="numeric",
                              default=101,
                              action="store",
                              dest="window_length",
                              metavar="Window Length",
                              help=paste("Length of the window for the moving average ",
                                         "(smoothing). Should be an odd integer.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--smooth_method"),
                              type="character",
                              default="pyramidinal",
                              action="store",
                              dest="smooth_method",
                              metavar="Smoothing Method",
                              help=paste("Method to use for smoothing: c(runmeans,pyramidinal)",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--num_ref_groups"),
                              type="numeric",
                              default=NULL,
                              action="store",
                              dest="num_ref_groups",
                              metavar="Number of reference groups",
                              help=paste("The number of reference groups or a list of",
                                         "indices for each group of reference indices in",
                                         "relation to reference_obs. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--ref_subtract_use_mean_bounds"),
                              type="logical",
                              default=TRUE,
                              action="store_false",
                              dest="ref_subtract_use_mean_bounds",
                              metavar="Reference Subtract use Mean Bounds",
                              help=paste("Determine means separately for each ref group, ",
                                         "then remove intensities within bounds of means",
                                         "[Default %default]",
                                         "Otherwise, uses mean of the means across groups."))

pargs <- optparse::add_option(pargs, c("--cluster_by_groups"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="cluster_by_groups",
                              metavar="Cluster by Groups",
                              help=paste("If observations are defined according to groups ",
                                         "(ie. patients), each group of cells will be ",
                                         "clustered separately. ([Default %default]",
                                         ", instead will use k_obs_groups setting)"))

pargs <- optparse::add_option(pargs, c("--k_obs_groups"),
                              type="numeric",
                              default=1,
                              action="store",
                              dest="k_obs_groups",
                              metavar="K number of Observation groups",
                              help=paste("Number of groups in which to break the observations.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--hclust_method"),
                              type="character",
                              default="ward.D2",
                              action="store",
                              dest="hclust_method",
                              metavar="Hierarchical Clustering Method",
                              help=paste("Method used for hierarchical clustering of cells. ",
                                         "Valid choices are: \"ward.D\", \"ward.D2\", \"single\"",
                                         ", \"complete\", \"average\", \"mcquitty\", \"median\", \"centroid\". ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--analysis_mode"),
                              type="character",
                              default="samples",
                              action="store",
                              dest="analysis_mode",
                              metavar="Analysis Mode",
                              help=paste("options(samples|subclusters|cells), ",
                                         "Grouping level for image filtering or HMM predictions.",
                                         "[Default %default] (fastest, but subclusters is ideal)"))

pargs <- optparse::add_option(pargs, c("--max_centered_threshold"),
                              type="numeric",
                              default=3,
                              action="store",
                              dest="max_centered_threshold",
                              metavar="Max Centered Threshold",
                              help=paste("The maximum value a value can have after",
                                         "centering. Also sets a lower bound of -1 * this value. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--scale_data"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="scale_data",
                              metavar="Scale Data",
                              help=paste("perform Z-scaling of logtransformed data ",
                                         "[Default %default]. ",
                                         "This may be turned on if you have very different ",
                                         "kinds of data for your normal and tumor samples. ",
                                         "For example, you need to use GTEx representative ", 
                                         "normal expression profiles rather than being able ",
                                         "to leverage normal single cell data that ",
                                         "goes with your experiment."))

pargs <- optparse::add_option(pargs, c("--HMM"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="HMM",
                              metavar="HMM",
                              help=paste("when set to True, runs HMM to predict CNV level. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--HMM_transition_prob"),
                              type="numeric",
                              default=1e-6,
                              action="store",
                              dest="HMM_transition_prob",
                              metavar="HMM Transition Probabiltie",
                              help=paste("transition probability in HMM",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--tumor_subcluster_pval"),
                              type="numeric",
                              default=0.01,
                              action="store",
                              dest="tumor_subcluster_pval",
                              metavar="Tumor Subcluster p-value",
                              help=paste("Max p-value for defining a significant tumor subcluster. ",
                                         "[Default %default]"))


pargs <- optparse::add_option(pargs, c("--tumor_subcluster_partition_method"),
                              type="character",
                              default="leiden",
                              action="store",
                              dest="tumor_subcluster_partition_method",
                              metavar="Tumor Subcluster Partition Method",
                              help=paste("c('leiden', 'random_trees', 'qnorm')",
                                         "[Default %default]",
                                         "method for defining tumor subclusters.",
                                         "leiden: Runs a nearest neighbor search, where communities are then partitionned with the Leiden algorithm.",
                                         "random_trees: Slow, uses permutation statistics w/ tree construction.",
                                         "qnorm: defines tree height based on the quantile defined by the tumor_subcluster_pval"))

pargs <- optparse::add_option(pargs, c("--HMM_report_by"),
                              type="character",
                              default="subcluster",
                              action="store",
                              dest="HMM_report_by",
                              metavar="HMM report by",
                              help=paste("c(cell, consensus, subcluster)",
                                         "[Default %default]",
                                         "Note, reporting is performed entirely",
                                         " separately from the HMM prediction.  ",
                                         "So, you can predict on subclusters, but ",
                                         "get per-cell level reporting (more voluminous output)."))

pargs <- optparse::add_option(pargs, c("--HMM_type"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="HMM_type",
                              metavar="HMM type",
                              help=paste("HMM model type. Options: (i6 or i3):",
                                         "i6: infercnv 6-state model (0, 0.5, 1,",
                                         " 1.5, 2, >2) where state emissions are ",
                                         "calibrated based on simulated CNV levels.\n",
                                         "i3: infercnv 3-state model (del, neutral, amp) ",
                                         "configured based on normal cells and HMM_i3_z_pval.\n",
                                         "[Default %default]"))

# pargs <- optparse::add_option(pargs, c("--HMM_i3_z_pval"),
#                               type="numeric",
#                               default=0.05,
#                               action="store",
#                               dest="HMM_i3_z_pval",
#                               metavar="HMM i3 z p-value",
#                               help=paste("p-value for HMM i3 state overlap",
#                                          "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--denoise"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="denoise",
                              metavar="Denoise",
                              help=paste("If True, turns on denoising according to options below",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--noise_filter"),
                              type="numeric",
                              default=NA,
                              action="store",
                              dest="noise_filter",
                              metavar="Noise Filter",
                              help=paste("Values +- from the reference cell mean will ",
                                         "be set to zero (whitening effect)",
                                         "[Default %default, instead will use ",
                                         "sd_amplifier below.]"))

pargs <- optparse::add_option(pargs, c("--sd_amplifier"),
                              type="numeric",
                              default=1.0,
                              action="store",
                              dest="sd_amplifier",
                              metavar="SD denoise amplifier",
                              help=paste("Noise is defined as mean(reference_cells) ",
                                         "+- sdev(reference_cells) * sd_amplifier ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--noise_logistic"),
                              type="logical",
                              default=TRUE,
                              action="store_false",
                              dest="noise_logistic",
                              metavar="Noise Logistic",
                              help=paste("use the noise_filter or sd_amplifier ",
                                         "based threshold (whichever is invoked) ",
                                         "as the midpoint in alogistic model for ",
                                         "downscaling values close to the mean. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--outlier_method_bound"),
                              type="character",
                              default="average_bound",
                              action="store",
                              dest="outlier_method_bound",
                              metavar="Outlier Method Bound",
                              help=paste("Method to use for bounding outlier values. ",
                                         "[Default %default]",
                                         "Will preferentially use outlier_lower_bound ",
                                         "and outlier_upper_bound if set."))

pargs <- optparse::add_option(pargs, c("--outlier_lower_bound"),
                              type="numeric",
                              default=NA,
                              action="store",
                              dest="outlier_lower_bound",
                              metavar="Outlier Lower Bound",
                              help=paste("Outliers below this lower bound ",
                                         "will be set to this value.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--outlier_upper_bound"),
                              type="numeric",
                              default=NA,
                              action="store",
                              dest="outlier_upper_bound",
                              metavar="Outlier Upper Bound",
                              help=paste("Outliers above this upper bound ",
                                         "will be set to this value.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--plot_steps"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="plot_steps",
                              metavar="Plot Steps",
                              help=paste("If true, saves infercnv objects and ",
                                         "plots data at the intermediate steps.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--final_scale_limits"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="final_scale_limits",
                              metavar="Final Scale Limits",
                              help=paste("The scale limits for the final heatmap ",
                                         "output by the run() method. ",
                                         "[Default %default] ",
                                         " Alt, c(low,high)"))

pargs <- optparse::add_option(pargs, c("--final_center_val"),
                              type="numeric",
                              default=NULL,
                              action="store",
                              dest="final_center_val",
                              metavar="Final Center Value",
                              help=paste("Center value for final heatmap output ",
                                         "by the run() method.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--debug"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="debug",
                              metavar="Debug",
                              help=paste("If true, output debug level logging.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--num_threads"),
                              type="numeric",
                              default=4,
                              action="store",
                              dest="num_threads",
                              metavar="Number of Threads",
                              help=paste("(int) number of threads for parallel steps. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--raw_counts_matrix"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="raw_counts_matrix",
                              metavar="Raw Counts Expression Data",
                              help=paste("the matrix of genes (rows) vs. cells (columns) ",
                                         "containing the raw counts. It'll be read via read.table()",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--gene_order_file"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="gene_order_file",
                              metavar="Gene Order File",
                              help=paste("data file containing the positions of ",
                                         "each gene along each chromosome in the genome. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--annotations_file"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="annotations_file",
                              metavar="Annotation File",
                              help=paste("a description of the cells, indicating ",
                                         "the cell type classifications. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--ref_group_names"),
                              type="character",
                              default=NULL,
                              action="store",
                              dest="ref_group_names",
                              metavar="Reference Groups Names",
                              help=paste("Names of groups from raw_counts_matrix whose cells",
                                         "are to be used as reference groups.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--delim"),
                              type="character",
                              action="store",
                              default="\t",
                              dest="delim",
                              metavar="Delimiter",
                              help=paste("Delimiter for reading expression matrix",
                                         " and writing matrices output.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--max_cells_per_group"),
                              type="numeric",
                              default=NULL,
                              action="store",
                              dest="max_cells_per_group",
                              metavar="Max Cells per group",
                              help=paste("maximun number of cells to use per group. ",
                                         "[Default %default] using all cells defined ",
                                         "in the annotations_file. This option is useful ",
                                         "for randomly subsetting the existing data ",
                                         "for a quicker preview run, such as using ",
                                         "50 cells per group instead of hundreds."))

pargs <- optparse::add_option(pargs, c("--log_file"),
                              type="character",
                              action="store",
                              default=NA,
                              dest="log_file",
                              metavar="Log",
                              help=paste("File for logging. If not given,",
                                         "logging will occur to console.",
                                         "[Default %default]"))

#pargs <- optparse::add_option(pargs, c("--contig_lab_size"),
#                              type="integer",
#                              action="store",
#                              default=1,
#                              dest="contig_label_size",
#                              metavar="Contig_Label_Size",
#                              help=paste("Used to increase or decrease the text labels",
#                                         "for the X axis (contig names).",
#                                         "[Default %default]"))

#pargs <- optparse::add_option(pargs, c("--color_safe"),
#                              type="logical",
#                              default=FALSE,
#                              action="store_true",
#                              dest="use_color_safe",
#                              metavar="Color_Safe",
#                              help=paste("To support the needs of those who see ",
#                                         "colors differently, use this option to",
#                                         "change the colors to a palette visibly ",
#                                         "distinct to all color blindness. ",
#                                         " [Default %default]"))

#pargs <- optparse::add_option(pargs, c("--title"),
#                              type="character",
#                              default="Copy Number Variation Inference",
#                              action="store",
#                              dest="fig_main",
#                              metavar="Figure_Title",
#                              help=paste("Title of the figure.",
#                                         "[Default %default]"))

#pargs <- optparse::add_option(pargs, c("--title_obs"),
#                              type="character",
#                              default="Observations (Cells)",
#                              action="store",
#                              dest="obs_main",
#                              metavar="Observations_Title",
#                              help=paste("Title of the observations matrix Y-axis.",
#                                         "[Default %default]"))

#pargs <- optparse::add_option(pargs, c("--title_ref"),
#                              type="character",
#                              default="References (Cells)",
#                              action="store",
#                              dest="ref_main",
#                              metavar="References_Title",
#                              help=paste("Title of the references matrix Y-axis (if used).",
#                                         "[Default %default]"))

#pargs <- optparse::add_option(pargs, c("--ngchm"),
#                              type="logical",
#                              action="store_true",
#                              default=FALSE,
#                              dest="ngchm",
#                              metavar="NextGen_HeatMap",
#                              help=paste("Create a Next Generation Clustered Heat Map"))

#pargs <- optparse::add_option(pargs, c("--path_to_shaidyMapGen"),
#                              type="character",
#                              action="store",
#                              default=NULL,
#                              dest="path_to_shaidyMapGen",
#                              metavar="Path_To_ShaidyMapGenp",
#                              help=paste("This is the pathway to the java application ShaidyMapGen.jar.",
#                                    "If this is not assigned, then an enviornmental variable that ",
#                                    "contains the "))

#pargs <- optparse::add_option(pargs, c("--gene_symbol"),
#                              type="character",
#                              action="store",
#                              default=NULL,
#                              dest="gene_symbol",
#                              metavar="Gene_Symbol",
#                              help=paste("The labeling type used to represent the genes in the expression",
#                                   "data. This needs to be passed in order to add linkouts to the ",
#                                   "genes. Possible gene label types to choose from are specified on",
#                                   "the broadinstitute/inferCNV wiki and bmbroom/NGCHM-config-biobase."))

pargs <- optparse::add_option(pargs, c("--no_plot"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="no_plot",
                              metavar="No Plot",
                              help=paste("don't make any of the images.",
                                         "Instead, generate all non-image outputs as part of the run.",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--no_prelim_plot"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="no_prelim_plot",
                              metavar="No Preliminary Plot",
                              help=paste("don't make the preliminary infercnv image",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--output_format"),
                              type="character",
                              default="png",
                              action="store",
                              dest="output_format",
                              metavar="Output Format",
                              help=paste("Output format for the figure. Default is NA, ",
                                          "which means to only write the text outputs ",
                                          "without generating the figure itself. ",
                                          "Other choices are \"png\" and \"pdf\"."))

pargs <- optparse::add_option(pargs, c("--median_filter"),
                              type="logical",
                              default=FALSE,
                              action="store_true",
                              dest="median_filter",
                              metavar="Median Filter",
                              help=paste("If True, turns on additional median",
                                         " filtering for an additional plot. ",
                                         "[Default %default]"))

pargs <- optparse::add_option(pargs, c("--top_n"),
                              type="numeric",
                              default=10,
                              action="store",
                              dest="top_n",
                              metavar="Top n exported CNV",
                              help=paste("(int) number of top CNVs to export with add_to_seurat.",
                                         "[Default %default]"))

args <- optparse::parse_args(pargs)

# Check arguments
#args <- check_arguments(args)

if (!is.null(args$final_scale_limits)) {
    if (grepl(',', args$final_scale_limits)) {
        args$final_scale_limits = as.double(strsplit(args$final_scale_limits, ","))
    }
}

if (!is.null(args$ref_group_names)) {
    args$ref_group_names = strsplit(args$ref_group_names, ",")[[1]]
} else {
    args$ref_group_names = c()
}


# Parse bounds
#bounds_viz <- c(NA,NA)
#if (!is.na(args$bound_threshold_vis)){
#    bounds_viz <- as.numeric(unlist(strsplit(args$bound_threshold_vis,",")))
#}
#if (length(bounds_viz) != 2){
#    error_message <- paste("Please use the correct format for the argument",
#                           "--vis_bound_threshold . Two numbers seperated",
#                           "by a comma is expected (lowerbound,upperbound)",
#                           ". As an example, to indicate that outliers are",
#                           "outside of -1 and 1 give the following.",
#                           "--vis_bound_threshold -1,1")
#    stop(error_message)
#}

# Set up logging file
#logging::basicConfig(level=args$log_level)
#if (!is.na(args$log_file)){
#    logging::addHandler(logging::writeToFile,
#                        file=args$log_file,
#                        level=args$log_level)
#}

# Log the input parameters
flog.info(paste("::Input arguments. Start."))
for (arg_name in names(args)){
    flog.info(paste(":Input_Argument:",arg_name,"=",args[[arg_name]],
                           sep=""))
}
flog.info(paste("::Input arguments. End."))

infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=args$raw_counts_matrix,
                                               gene_order_file=args$gene_order_file,
                                               annotations_file=args$annotations_file,
                                               ref_group_names=args$ref_group_names,
                                               delim=args$delim,
                                               max_cells_per_group=args$max_cells_per_group,
                                               chr_exclude=args$chr_exclude)

infercnv_obj = infercnv::run(infercnv_obj=infercnv_obj,
                            cutoff=args$cutoff,
                            min_cells_per_gene=args$min_cells_per_gene,
                            out_dir=args$out_dir,
                            analysis_mode=args$analysis_mode,
                            window_length=args$window_length,
                            smooth_method=args$smooth_method,
                            num_ref_groups=args$num_ref_groups,
                            ref_subtract_use_mean_bounds=args$ref_subtract_use_mean_bounds,
                            max_centered_threshold=args$max_centered_threshold,
                            tumor_subcluster_pval=args$tumor_subcluster_pval,
                            tumor_subcluster_partition_method=args$tumor_subcluster_partition_method,
                            HMM=args$HMM,
                            HMM_transition_prob=args$HMM_transition_prob,
                            HMM_report_by=args$HMM_report_by,
                            HMM_type=args$HMM_type,
                            # HMM_i3_z_pval=args$HMM_i3_z_pval,
                            #sim_method=args$sim_method,
                            #sim_foreground=args$sim_foreground,
                            scale_data=args$scale_data,
                            denoise=args$denoise,
                            noise_filter=args$noise_filter,
                            sd_amplifier=args$sd_amplifier,
                            noise_logistic=args$noise_logistic,
                            cluster_by_groups=args$cluster_by_groups,
                            k_obs_groups=args$k_obs_groups,
                            outlier_method_bound=args$outlier_method_bound,
                            outlier_lower_bound=args$outlier_lower_bound,
                            outlier_upper_bound=args$outlier_upper_bound,
                            hclust_method=args$hclust_method,
                            #remove_genes_at_chr_ends=args$remove_genes_at_chr_ends,
                            #mask_nonDE_genes=args$mask_nonDE_genes,
                            #mask_nonDE_pval=args$mask_nonDE_pval,
                            #test.use=args$test.use,
                            #require_DE_all_normals=args$require_DE_all_normals,
                            plot_steps=args$plot_steps,
                            no_plot=args$no_plot,
                            no_prelim_plot=args$no_prelim_plot,
                            output_format=args$output_format,
                            debug=args$debug,
                            #prune_outliers=args$prune_outliers,
                            final_scale_limits=args$final_scale_limits,
                            final_center_val=args$final_center_val,
                            #reuse_subtracted=args$reuse_subtracted,
                            num_threads=args$num_threads#,
                            #hspike_aggregate_normals =args$hspike_aggregate_normals
                            )

if (args$median_filter) {

    infercnv_obj = infercnv::apply_median_filtering(infercnv_obj)
    
    if (is.null(args$final_scale_limits)) {
        args$final_scale_limits = "auto"
    }
    if (is.null(args$final_center_val)) {
        args$final_center_val = 1
    }
    
    plot_cnv(infercnv_obj,
             k_obs_groups=args$k_obs_groups,
             cluster_by_groups=args$cluster_by_groups,
             out_dir=args$out_dir,
             x.center=args$final_center_val,
             x.range=args$final_scale_limits,
             title="inferCNV",
             output_filename="infercnv_pdf",
             write_expr_matrix=TRUE)

}

if (args$HMM) {
    add_to_seurat(infercnv_output_path=args$out_dir,
                  top_n=args$top_n
    )
}
