#' MCMC infercnv_combined class
#' 
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of three/six Copy Number Variation states 
#' (i6 states: 1-2,Deletion; 3,Neutral; 4-6,Amplification; 7,cnLOH,if enable_cnLOH),
#' (i3 states: 1,Deletion; 2,Neutral; 3,Amplification; 4,cnLOH,if enable_cnLOH),
#' for CNV's identified by inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster 
#' and each individual cell line in the CNV.
#' 
#' @slot bugs_model BUGS model.
#' 
#' @slot cell_gene List containing the Cells and Genes that make up each CNV.
#' 
#' @slot cnv_probabilities Probabilities of each CNV belonging to a particular state from 0 (least likely)to 1 (most likely).
#' 
#' @slot cell_probabilities Probabilities of each cell being in a particular state, from 0 (least likely)to 1 (most likely).
#' 
#' @slot cnv_regions ID for each CNV found by the HMM.
#' 
#' @slot infercnv_obj InferCNV obj.
#' 
#' @slot infercnv_obj_mu Mean values to be used for determining the distribution of each cell line.
#' 
#' @slot infercnv_obj_sig fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line.
#' 
#' @slot infercnv_allele_obj InferCNV_allele obj.
#' 
#' @export
#' 
MCMC_infercnv_combined <- methods::setClass("MCMC_infercnv_combined", 
                                            slots = c(bugs_model = "character",
                                                      cell_gene = "list",
                                                      cnv_probabilities = "list",
                                                      cell_probabilities = "list",
                                                      cnv_regions = "factor",
                                                      # infercnv_obj = "infercnv",
                                                      infercnv_obj_mu = "numeric",
                                                      infercnv_obj_sig = "numeric"),
                                            contains = "infercnv_allele")

# file_path the path to HMM reports
# file_token file name
# infercnv based obj
# infercnv_allele based obj
# two modes: snp and gene
initialize_combined_mcmc <- function(infercnv_allele_obj,
                                     file_path,
                                     file_token,
                                     # infercnv_obj,
                                     infercnv_obj_mu,
                                     infercnv_obj_sig,
                                     mode = c("snp_level", "gene_level"),
                                     HMM_type = c("i6", "i3"),
                                     enable_cnLOH = TRUE) {
  
    mode <- match.arg(mode)
    HMM_type <- match.arg(HMM_type)
    
    ## Load the files for cnv predictions
    cell_groups_df <- read.table(file.path(file_path, paste0(file_token,".cell_groupings")), 
                                 header = T,
                                 check.names = FALSE,
                                 sep = "\t")
    pred_cnv_genes_df <- read.table(file.path(file_path,
        paste0(file_token,".pred_cnv_genes.dat")),
                                    header = T,
                                    check.names = FALSE,
                                    sep = "\t",
                                    stringsAsFactors = TRUE)
    
    flog.info(sprintf("Initializing MCMC_infercnv_combined obj at %s ...", mode))
    
    MCMC_inferCNV_combined_obj <- new("MCMC_infercnv_combined", infercnv_allele_obj)

    MCMC_inferCNV_combined_obj@infercnv_obj_mu = infercnv_obj_mu
    MCMC_inferCNV_combined_obj@infercnv_obj_sig = infercnv_obj_sig
    MCMC_inferCNV_combined_obj@infercnv_obj_mu = infercnv_obj_mu

    # mcmc_combined <- new("MCMC_infercnv_combined",
    #                       infercnv_obj = infercnv_obj,
    #                       infercnv_obj_mu = infercnv_obj_mu,
    #                       infercnv_obj_sig = infercnv_obj_sig,
    #                       infercnv_allele_obj = infercnv_allele_obj)
    
    MCMC_inferCNV_combined_obj@bugs_model <- ifelse(HMM_type == "i6",
                                                    ifelse(mode == "snp_level",
                                                           ifelse(enable_cnLOH,
                                                                  system.file("BUGS_Combined_Model_i7", package = "infercnv"),
                                                                  system.file("BUGS_Combined_Model_i6", package = "infercnv")),
                                                           ifelse(enable_cnLOH,
                                                                  system.file("BUGS_SNP2Gene_Combined_Model_i7", package = "infercnv"),
                                                                  system.file("BUGS_SNP2Gene_Combined_Model_i6", package = "infercnv"))),
                                                    ifelse(mode == "snp_level",
                                                           ifelse(enable_cnLOH,
                                                                  system.file("BUGS_Combined_Model_i4", package = "infercnv"),
                                                                  system.file("BUGS_Combined_Model_i3", package = "infercnv")),
                                                           ifelse(enable_cnLOH,
                                                                  system.file("BUGS_SNP2Gene_Combined_Model_i4", package = "infercnv"),
                                                                  system.file("BUGS_SNP2Gene_Combined_Model_i3", package = "infercnv"))))

    MCMC_inferCNV_combined_obj@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
    
    MCMC_inferCNV_combined_obj <- getGenesCells_combined(MCMC_inferCNV_combined_obj,
                                                         infercnv_allele_obj,
                                                         pred_cnv_genes_df,
                                                         cell_groups_df,
                                                         mode = mode)

    return(MCMC_inferCNV_combined_obj)

    # if(mode == "snp_level"){
    #   
    #   flog.info("initializing MCMC_infercnv_combined obj at snp level ...")
    #   
    #   mcmc_combined_snp <- new("MCMC_infercnv_combined",
    #                            infercnv_obj = infercnv_obj,
    #                            infercnv_obj_mu = infercnv_obj_mu,
    #                            infercnv_obj_sig = infercnv_obj_sig,
    #                            infercnv_allele_obj = infercnv_allele_obj)
    #   mcmc_combined_snp@bugs_model <- ifelse(HMM_type == "i6",
    #                                          system.file("BUGS_Combined_Model_i6",package = "infercnv"),
    #                                          system.file("BUGS_Combined_Model_i3",package = "infercnv"))
    #     
    #   mcmc_combined_snp@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
    #   
    #   mcmc_combined_snp <- getGenesCells_combined(mcmc_combined_snp,
    #                                               pred_cnv_genes_df, 
    #                                               cell_groups_df,
    #                                               mode = mode)
    #   return(mcmc_combined_snp)
    #   
    # } else if(mode == "gene_level"){
    #   flog.info("initializing MCMC_infercnv_combined obj at gene level ...")
    #   
    #   mcmc_combined_gene <- new("MCMC_infercnv_combined",
    #                             infercnv_obj = infercnv_obj,
    #                             infercnv_obj_mu = infercnv_obj_mu,
    #                             infercnv_obj_sig = infercnv_obj_sig,
    #                             infercnv_allele_obj = infercnv_allele_obj)
    #   mcmc_combined_gene@bugs_model <- ifelse(HMM_type == "i6",
    #                                           system.file("BUGS_SNP2Gene_Combined_Model_i6",package = "infercnv"),
    #                                           system.file("BUGS_SNP2Gene_Combined_Model_i3",package = "infercnv"))
    # 
    #   mcmc_combined_gene@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
    #   
    #   mcmc_combined_gene <- getGenesCells_combined(mcmc_combined_gene,
    #                                                pred_cnv_genes_df, 
    #                                                cell_groups_df,
    #                                                mode = mode)
    #   return(mcmc_combined_gene)
    # }
  
}

getGenesCells_combined <- function(MCMC_inferCNV_combined_obj,
                                   infercnv_allele_obj,
                                   pred_cnv_genes_df,
                                   cell_groups_df, 
                                   mode = c("snp_level", "gene_level")) {
    mode = match.arg(mode)
    
    MCMC_inferCNV_combined_obj@cell_gene <- lapply(MCMC_inferCNV_combined_obj@cnv_regions, function(x) {

        current_cnv <- pred_cnv_genes_df[which(x == pred_cnv_genes_df$gene_region_name),]
        genes <- current_cnv$gene
        
        infercnv_gene_idx <- which(row.names(infercnv_allele_obj@expr.data) %in% genes)
        
        # if(mode == "snp_level"){
        #   infercnv_allele_gene <- obj@infercnv_allele_obj@SNP_info$gene[obj@infercnv_allele_obj@SNP_info$gene %in% genes]
        #   infercnv_allele_gene_idx <- which(obj@infercnv_allele_obj@SNP_info$gene %in% 
        #                                       infercnv_allele_gene)
        # } else{
        #   infercnv_allele_gene_idx <- which(obj@infercnv_allele_obj@SNP_info$gene %in% genes)
        # }
        
        infercnv_allele_gene_idx <- which(MCMC_inferCNV_combined_obj@SNP_info$gene %in% genes)
        
        sub_cells <- unique(current_cnv$cell_group_name)
        infercnv_cells_idx <- which(colnames(infercnv_allele_obj@expr.data) %in% 
                                    cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
        infercnv_allele_cells_idx <- which(colnames(infercnv_allele_obj@allele.expr.data) %in% 
                                           cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
        state <- unique(current_cnv$state)
        
        if (length(infercnv_allele_gene_idx) == 0){
          
            if(mode == "snp_level"){
                return(list("cnv_regions" = x, 
                            
                            "infercnv_Genes" = infercnv_gene_idx,
                            "infercnv_Cells" = infercnv_cells_idx,
                            
                            "infercnv_allele_SNPs" = NULL, 
                            "infercnv_allele_Cells" = NULL, 
                            
                            "r.array" = NULL,
                            "n.sc.array" = NULL,
                            "l.array" = NULL,
                            "n.bulk.array" = NULL,
                            "I.j" = NULL,
                            
                            "State" = state))
            } else{
                return(list("cnv_regions" = x, 
                            
                            "infercnv_Genes" = infercnv_gene_idx,
                            "infercnv_Cells" = infercnv_cells_idx,
                            
                            "infercnv_allele_Genes" = NULL, 
                            "infercnv_allele_Cells" = NULL, 
                            
                            "State" = state))
            }
            
        }
        
        if (mode == "snp_level") {
          
            ## extract affected genes
            genes.of.interest <- unique(MCMC_inferCNV_combined_obj@SNP_info[infercnv_allele_gene_idx]$gene)
            
            ## associate each gene factor with a set of snps
            candidate_snp <- MCMC_inferCNV_combined_obj@SNP_info[infercnv_allele_gene_idx]
            genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
                names(candidate_snp)[which(candidate_snp$gene %in% genes.of.interest[i])]
            })
            names(genes2snps.dict) <- genes.of.interest
            
            multi_arrays <- set_array(r.maf = infercnv_allele_obj@allele.count.data[infercnv_allele_gene_idx,
                                                                                    infercnv_allele_cells_idx,
                                                                                    drop = F], 
                                      n.sc = infercnv_allele_obj@allele.coverage.data[infercnv_allele_gene_idx,
                                                                                      infercnv_allele_cells_idx,
                                                                                      drop = F],
                                      genes2snps.dict = genes2snps.dict, 
                                      numCells = length(infercnv_allele_cells_idx))
            
            return(list("cnv_regions" = x, 
                        
                        "infercnv_Genes" = infercnv_gene_idx,
                        "infercnv_Cells" = infercnv_cells_idx,
                        
                        "infercnv_allele_SNPs" = infercnv_allele_gene_idx, 
                        "infercnv_allele_Cells" = infercnv_allele_cells_idx, 
                        
                        
                        "r.array" = multi_arrays$r.array,
                        "n.sc.array" = multi_arrays$n.sc.array,
                        "l.array" = multi_arrays$l.array,
                        "n.bulk.array" = multi_arrays$n.bulk.array,
                        "I.j"=multi_arrays$I.j,
                        
                        "State" = state))
            
        } else {
            return(list("cnv_regions" = x, 
                        
                        "infercnv_Genes" = infercnv_gene_idx,
                        "infercnv_Cells" = infercnv_cells_idx,
                        
                        "infercnv_allele_Genes" = infercnv_allele_gene_idx, 
                        "infercnv_allele_Cells" = infercnv_allele_cells_idx, 
                        
                        "State" = state))
        }
    })
    return(MCMC_inferCNV_combined_obj)
}

#' @title inferCNVCombinedBayesNet: Run Bayesian Network Mixture Model Leveraging both
#' Expression and Allele Data To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of three/six Copy Number Variation states 
#' (i6 states: 1-2,Deletion; 3,Neutral; 4-6,Amplification; 7,cnLOH,if enable_cnLOH),
#' (i3 states: 1,Deletion; 2,Neutral; 3,Amplification; 4,cnLOH,if enable_cnLOH),
#' for CNV's identified by inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster 
#' and each individual cell line in the CNV.
#' 
#' @param combined_file_path Location of the directory of the combined HMM outputs.
#' 
#' @param combined_file_token (string) String token that contains some info on settings used to name combined-based files.
#' 
#' @param infercnv_obj InferCNV object.
#' 
#' @param infercnv_allele_obj InferCNV allele object.
#' 
#' @param allele_mode The type of allele data provided. "snp_level" or "gene_level".
#' 
#' @param method The way integrating HMM boundaries from two methods: union/common.
#' "union" keeps those expression-specific gene that are not found in allele data during modeling,
#' while "common" removes those ones.
#' 
#' @param HMM_type The type of combined model running in Bayesian Method. i6 or i3.
#' 
#' @param enable_cnLOH Option for detecting cnLOH events by adding one more state in Bayesian Method. (Default: True)
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param cores Option to run parallel by specifying the number of cores to be used. (Default: 5)
#' 
#' @return Returns a MCMC_infercnv_combined obj and posterior probability of being in one of three/six Copy Number Variation states
#' (i6 states: 1-2,Deletion; 3,Neutral; 4-6,Amplification; 7,cnLOH,if enable_cnLOH),
#' (i3 states: 1,Deletion; 2,Neutral; 3,Amplification; 4,cnLOH,if enable_cnLOH), for CNV's identified by inferCNV's HMM.
#'
#' @export
#' 
#' @examples 
#' data(infercnv_data_example)
#' data(infercnv_annots_example)
#' data(infercnv_genes_example)
#' data(infercnv_allele_alt_example)
#' data(infercnv_allele_tot_example)
#' data(infercnv_object_allele_gene_example)
#' 
#' infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix = infercnv_data_example,
#'                                                           raw_allele_matrix = infercnv_allele_alt_example,
#'                                                           raw_coverage_matrix = infercnv_allele_tot_example,
#'                                                           gene_order_file = infercnv_genes_example,
#'                                                           annotations_file = infercnv_annots_example,
#'                                                           ref_group_names=c("Microglia/Macrophage"),
#'                                                           snp_split_by = "::")
#'           
#' out_dir = tempfile()
#' file_gene_token <- "HMM_gene_pred"
#' infercnv_object_example <- infercnv::run(infercnv_object_example,
#'                                          cutoff=1,
#'                                          out_dir=out_dir, 
#'                                          cluster_by_groups=F, 
#'                                          denoise=FALSE,
#'                                          HMM=TRUE,
#'                                          analysis_mode = "samples",
#'                                          no_plot=TRUE)
#'                                          
#' hmm_allele_gene_obj_HMM_samples <- infercnv:::allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_object_allele_gene_example,
#'                                                                                                     trim = 0)
#' infercnv:::generate_cnv_region_reports(hmm_allele_gene_obj_HMM_samples, 
#'                                        output_filename_prefix=file_gene_token,
#'                                        out_dir=out_dir,
#'                                        ignore_neutral_state = 2,
#'                                        by="subcluster")
#'                                        
#' combined_obj <- infercnv::fusion_HMM_report_two_methods(expression_file_path = out_dir,
#'                                                         expression_file_token = "17_HMM_predHMMi6.hmm_mode-samples",
#'                                                         allele_file_path = out_dir,
#'                                                         allele_file_token = file_gene_token,
#'                                                         infercnv_obj = infercnv_object_example,
#'                                                         infercnv_allele_obj = infercnv_object_allele_gene_example,
#'                                                         method = "common",
#'                                                         HMM_type = "i6",
#'                                                         cluster_by_groups=F,
#'                                                         output_path = out_dir,
#'                                                         output_prefix = "combined")
#' mcmc_combined <- infercnv::inferCNVCombinedBayesNet(combined_file_path = out_dir,
#'                                                     combined_file_token = "combined",
#'                                                     infercnv_obj = infercnv_object_example,
#'                                                     infercnv_allele_obj = infercnv_object_allele_gene_example,
#'                                                     allele_mode = "gene_level",
#'                                                     method = "common",
#'                                                     HMM_type = "i6",
#'                                                     enable_cnLOH = TRUE,
#'                                                     output_path = file.path(out_dir,"bayesian_folder"))

inferCNVCombinedBayesNet <- function(combined_file_path,
                                     combined_file_token,
                                     # infercnv_obj,
                                     infercnv_allele_obj,
                                     allele_mode = c("snp_level", "gene_level"),
                                     method = c("common", "union"),
                                     HMM_type = c("i6", "i3"),
                                     enable_cnLOH = TRUE,
                                     output_path,
                                     cores = 5) {
  
    allele_mode <- match.arg(allele_mode)
    method <- match.arg(method)
    HMM_type <- match.arg(HMM_type)
    
    if (!dir.exists(file.path(output_path))) {
        dir.create(file.path(output_path), recursive = T)
        flog.info(paste("Creating the following Directory:", output_path))
    }
    
    if (HMM_type == "i6") {
        cnv_mean_sd <- get_spike_dists(infercnv_allele_obj@.hspike)
        mean <- c(cnv_mean_sd[["cnv:0.01"]]$mean,
                  cnv_mean_sd[["cnv:0.5"]]$mean,
                  cnv_mean_sd[["cnv:1"]]$mean,
                  cnv_mean_sd[["cnv:1.5"]]$mean,
                  cnv_mean_sd[["cnv:2"]]$mean,
                  cnv_mean_sd[["cnv:3"]]$mean)
        
        sd <- c(1/(cnv_mean_sd[["cnv:0.01"]]$sd^2),
                1/(cnv_mean_sd[["cnv:0.5"]]$sd^2),
                1/(cnv_mean_sd[["cnv:1"]]$sd^2),
                1/(cnv_mean_sd[["cnv:1.5"]]$sd^2),
                1/(cnv_mean_sd[["cnv:2"]]$sd^2),
                1/(cnv_mean_sd[["cnv:3"]]$sd^2))
    } else {
        suppressMessages(invisible(capture.output(cnv_mean_sd <- .i3HMM_get_sd_trend_by_num_cells_fit(infercnv_allele_obj))))
        mean <- c(cnv_mean_sd$mu - cnv_mean_sd$mean_delta,
                  cnv_mean_sd$mu,
                  cnv_mean_sd$mu + cnv_mean_sd$mean_delta)
        
        sd <- c(1/(cnv_mean_sd$sigma^2),
                1/(cnv_mean_sd$sigma^2),
                1/(cnv_mean_sd$sigma^2))
    }
    
    if (method == "common") {
        infercnv_allele_obj <- remove_genes(infercnv_obj,
                                            which(! rownames(infercnv_allele_obj@gene_order) %in% 
                                                    infercnv_allele_obj@SNP_info$gene))
    }
    
    MCMC_inferCNV_combined_obj <- initialize_combined_mcmc(infercnv_allele_obj = infercnv_allele_obj,
                                                           file_path = combined_file_path,
                                                           file_token = combined_file_token,
                                                           # infercnv_obj = infercnv_obj,
                                                           infercnv_obj_mu = mean,
                                                           infercnv_obj_sig = sd,
                                                           mode = allele_mode,
                                                           HMM_type = HMM_type,
                                                           enable_cnLOH = enable_cnLOH)
                                                         
    
    flog.info(paste("The number of affected regions:", length(MCMC_inferCNV_combined_obj@cell_gene)))
    
    flog.info(sprintf("Start running Gibbs sampling in %s mode", HMM_type))
    
    flog.info(sprintf("Enable identifing cnLOH event: %s", enable_cnLOH))
    
    start_time <- Sys.time()
    
    mcmc <- mclapply(MCMC_inferCNV_combined_obj@cell_gene, function(x) {

        if (allele_mode == "snp_level") {
            # if(is.null(x$r.array)){
            #   MCMC_inferCNV_combined_obj@bugs_model <- ifelse(type == "i6",
            #                                                   system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
            #                                                   system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv"))
            # } 
            samples <- run_combined_snp_mcmc(bugs = ifelse(is.null(x$r.array),
                                                           ifelse(HMM_type == "i6",
                                                                  system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
                                                                  system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv")),
                                                           MCMC_inferCNV_combined_obj@bugs_model),
                                             HMM_type = HMM_type,
                                             enable_cnLOH = ifelse(is.null(x$r.array),
                                                                   FALSE,
                                                                   enable_cnLOH),
                                             r.array = x$r.array,
                                             n.sc.array = x$n.sc.array,
                                             l.array = x$l.array,
                                             n.bulk.array = x$n.bulk.array,
                                             nAlleleGenes = length(x$I.j),
                                             nGenes = length(x$infercnv_Genes),
                                             nCells = length(x$infercnv_Cells),
                                             I.j = x$I.j,
                                             gexp = MCMC_inferCNV_combined_obj@expr.data[x$infercnv_Genes,
                                                                                         x$infercnv_Cells,
                                                                                         drop = F],
                                             mu = MCMC_inferCNV_combined_obj@infercnv_obj_mu,
                                             sig = MCMC_inferCNV_combined_obj@infercnv_obj_sig,
                                             pe = 0.1,
                                             mono = 0.7)
            
            # plot_mcmc(samples = samples,
            #           region = x$cnv_regions,
            #           output_path = output_path)
          
        } else {
          
            # if(is.null(x$infercnv_allele_Genes)){
            #   MCMC_inferCNV_combined_obj@bugs_model <- ifelse(type == "i6",
            #                                                   system.file("BUGS_Mixture_Model_i6_gene_test",package = "infercnv"),
            #                                                   system.file("BUGS_Mixture_Model_i3_gene_test",package = "infercnv"))
            # } 
            samples <- run_combined_gene_mcmc(bugs = ifelse(is.null(x$infercnv_allele_Genes),
                                                            ifelse(HMM_type == "i6",
                                                                   system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
                                                                   system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv")),
                                                            MCMC_inferCNV_combined_obj@bugs_model),
                                              HMM_type = HMM_type,
                                              enable_cnLOH = ifelse(is.null(x$infercnv_allele_Genes),
                                                                    FALSE,
                                                                    enable_cnLOH),
                                              r = MCMC_inferCNV_combined_obj@allele.count.data[x$infercnv_allele_Genes,
                                                                                               x$infercnv_allele_Cells,
                                                                                               drop = F],
                                              n.sc = MCMC_inferCNV_combined_obj@allele.coverage.data[x$infercnv_allele_Genes,
                                                                                                     x$infercnv_allele_Cells,
                                                                                                     drop = F],
                                              nAlleleGenes = length(x$infercnv_allele_Genes),
                                              nGenes = length(x$infercnv_Genes),
                                              nCells = length(x$infercnv_Cells),
                                              gexp = MCMC_inferCNV_combined_obj@expr.data[x$infercnv_Genes,
                                                                                          x$infercnv_Cells,
                                                                                          drop = F],
                                              mu = MCMC_inferCNV_combined_obj@infercnv_obj_mu,
                                              sig = MCMC_inferCNV_combined_obj@infercnv_obj_sig,
                                              pe = 0.1,
                                              mono = 0.7)
            
            # plot_mcmc(samples = samples,
            #           region = x$cnv_regions,
            #           output_path = output_path)
        }
    }, mc.cores = cores)
    #})
    end_time <- Sys.time()
    flog.info(paste0("MCMC running time: ",
                     difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
    
    flog.info("Retrieving the posterior probabilities of CNVs ")
    MCMC_inferCNV_combined_obj <- get_posterior_prob(MCMC_inferCNV_obj = MCMC_inferCNV_combined_obj, mcmc = mcmc)
    
    saveRDS(MCMC_inferCNV_combined_obj, file = file.path(output_path, "MCMC_inferCNV_combined_obj.rds"))
    
    flog.info("Plotting the distribution of posterior probabilities")
    plot_posterior_prob(MCMC_inferCNV_obj = MCMC_inferCNV_combined_obj, output_path)
    
    flog.info("Plotting diagnostic statistics")
    plot_Diagnostics(mcmc = mcmc, output_path)
    
    return(MCMC_inferCNV_combined_obj)
}

run_combined_snp_mcmc <- function(bugs,
                                  HMM_type = c("i6", "i3"),
                                  enable_cnLOH,
                                  r.array,
                                  n.sc.array,
                                  l.array,
                                  n.bulk.array,
                                  nAlleleGenes,
                                  nGenes,
                                  nCells,
                                  I.j,
                                  gexp,
                                  mu,
                                  sig,
                                  pe = 0.1,
                                  mono = 0.7) {
  
    HMM_type <- match.arg(HMM_type)
    
    #browser()
    input <- list('r' = r.array,
                  'n.sc' = n.sc.array,
                  'l' = l.array,
                  'n.bulk' = n.bulk.array,
                  
                  'J' = nAlleleGenes,  # how many genes (snp)
                  'JJ' = nGenes, # how many genes (gene)
                  'K' = nCells, # how many cells
                  'I.j' = I.j,
                  
                  'gexp' = gexp,
                  'mu' = mu,
                  'sig' = sig,
                  
                  'pseudo' = pe,
                  'mono' = mono)
    
    if (HMM_type == "i6") {
        if (enable_cnLOH) {
            inits <- list(list(epsilon = rep(1, input$K)),
                          list(epsilon = rep(2, input$K)),
                          list(epsilon = rep(3, input$K)),
                          list(epsilon = rep(4, input$K)),
                          list(epsilon = rep(5, input$K)),
                          list(epsilon = rep(6, input$K)),
                          list(epsilon = rep(7, input$K)))
        } else {
            inits <- list(list(epsilon = rep(1, input$K)),
                          list(epsilon = rep(2, input$K)),
                          list(epsilon = rep(3, input$K)),
                          list(epsilon = rep(4, input$K)),
                          list(epsilon = rep(5, input$K)),
                          list(epsilon = rep(6, input$K)))
        }
    } else {
        if (enable_cnLOH){
            inits <- list(list(epsilon = rep(1, input$K)),
                          list(epsilon = rep(2, input$K)),
                          list(epsilon = rep(3, input$K)),
                          list(epsilon = rep(4, input$K)))
        } else {
            inits <- list(list(epsilon = rep(1, input$K)),
                          list(epsilon = rep(2, input$K)),
                          list(epsilon = rep(3, input$K)))
        }
    }
    
    
    model <- rjags::jags.model(bugs,
                               data = input,
                               inits = inits,
                               n.chains = length(inits), 
                               n.adapt = 500,
                               quiet = F)
    update(model, 200, quiet = F)
    parameters <- c('theta', 'epsilon')
    samples <- rjags::coda.samples(model,
                                   parameters, 
                                   n.iter=1e3, 
                                   thin = 10,
                                   quiet = F)
    
    return(samples)
}

run_combined_gene_mcmc <- function(bugs,
                                   HMM_type = c("i6","i3"),
                                   enable_cnLOH,
                                   r,
                                   n.sc,
                                   nAlleleGenes,
                                   nGenes,
                                   nCells,
                                   gexp,
                                   mu,
                                   sig,
                                   pe = 0.1,
                                   mono = 0.7) {
  
  HMM_type <- match.arg(HMM_type)
  
  #browser()
  input <- list('r' = r,
                'n.sc' = n.sc,
                'l' = rowSums(r > 0),
                'n.bulk' = rowSums(n.sc > 0),
                
                'J' = nAlleleGenes,  # how many genes (snp)
                'JJ' = nGenes, # how many genes (gene)
                'K' = nCells, # how many cells
                
                'gexp' = gexp,
                'mu' = mu,
                'sig' = sig,
                
                'pseudo' = pe,
                'mono' = mono)
  
  if(HMM_type == "i6"){
    if(enable_cnLOH){
      inits <- list(
        list(epsilon = rep(1, input$K)),
        list(epsilon = rep(2, input$K)),
        list(epsilon = rep(3, input$K)),
        list(epsilon = rep(4, input$K)),
        list(epsilon = rep(5, input$K)),
        list(epsilon = rep(6, input$K)),
        list(epsilon = rep(7, input$K)))
    } else{
      inits <- list(
        list(epsilon = rep(1, input$K)),
        list(epsilon = rep(2, input$K)),
        list(epsilon = rep(3, input$K)),
        list(epsilon = rep(4, input$K)),
        list(epsilon = rep(5, input$K)),
        list(epsilon = rep(6, input$K)))
    }
  } else{
    if(enable_cnLOH){
      inits <- list(
        list(epsilon = rep(1, input$K)),
        list(epsilon = rep(2, input$K)),
        list(epsilon = rep(3, input$K)),
        list(epsilon = rep(4, input$K)))
    } else{
      inits <- list(
        list(epsilon = rep(1, input$K)),
        list(epsilon = rep(2, input$K)),
        list(epsilon = rep(3, input$K)))
    }
  }
  
  model <- rjags::jags.model(bugs,
                             data = input,
                             inits = inits,
                             n.chains = length(inits), 
                             n.adapt = 500,
                             quiet = F)
  update(model, 200, quiet = F)
  parameters <- c('theta', 'epsilon')
  
  samples <- rjags::coda.samples(model,
                                 parameters, 
                                 n.iter = 1e3, 
                                 thin = 10,
                                 quiet = F)
  
  return(samples)
  
}

#' @title fusion_HMM_report_two_methods 
#'
#' @description This function aims to integrate HMM boundaries from expression/allele based methods.
#' 
#' @param expression_file_path Location of the directory of the inferCNV outputs.
#' 
#' @param expression_file_token (string) String token that contains some info on settings used to name expression-based files.
#' 
#' @param allele_file_path Location of the directory of the inferCNV_allele outputs.
#' 
#' @param allele_file_token (string) String token that contains some info on settings used to name allele-based files.
#' 
#' @param infercnv_obj InferCNV object.
#' 
#' @param infercnv_allele_obj InferCNV allele object.
#' 
#' @param method The way integrating HMM boundaries from two methods: union/common.
#' "union" keeps those expression-specific gene that are not found in allele data during modeling,
#' while "common" removes those ones.
#' 
#' @param HMM_type The type of combined model running in Bayesian Method. i6 or i3.
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param output_prefix name prefix when naming output files.
#' 
#' @param HMM_report_by cell, consensus, subcluster (default: subcluster)
#' 
#' @param ... Additional parameters passed to infercnv::plot_cnv()
#' 
#' @return Returns an Infercnv hmm-based obj, and HMM states store in the slot of expr.data
#' (i6 states: 1-2,Deletion; 3,Neutral; 4-6,Amplification; 7,cnLOH; 8, Opposed prediction),
#' (i3 states: 1,Deletion; 2,Neutral; 3,Amplification; 4,cnLOH; 5, Opposed prediction), for CNV's identified by inferCNV's HMM.
#'
#' @export
#' 
#' @examples 
#' data(infercnv_data_example)
#' data(infercnv_annots_example)
#' data(infercnv_genes_example)
#' data(infercnv_allele_alt_example)
#' data(infercnv_allele_tot_example)
#' data(infercnv_object_allele_gene_example)
#' 
#' infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix = infercnv_data_example,
#'                                                           raw_allele_matrix = infercnv_allele_alt_example,
#'                                                           raw_coverage_matrix = infercnv_allele_tot_example,
#'                                                           gene_order_file = infercnv_genes_example,
#'                                                           annotations_file = infercnv_annots_example,
#'                                                           ref_group_names=c("Microglia/Macrophage"),
#'                                                           snp_split_by = "::")
#'           
#' out_dir = tempfile()
#' file_gene_token <- "HMM_gene_pred"
#' infercnv_object_example <- infercnv::run(infercnv_object_example,
#'                                          cutoff=1,
#'                                          out_dir=out_dir, 
#'                                          cluster_by_groups=F, 
#'                                          denoise=FALSE,
#'                                          HMM=TRUE,
#'                                          analysis_mode = "samples",
#'                                          no_plot=TRUE)
#'                                          
#' hmm_allele_gene_obj_HMM_samples <- infercnv:::allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_object_allele_gene_example,
#'                                                                                                     trim = 0)
#' infercnv:::generate_cnv_region_reports(hmm_allele_gene_obj_HMM_samples, 
#'                                        output_filename_prefix=file_gene_token,
#'                                        out_dir=out_dir,
#'                                        ignore_neutral_state = 2,
#'                                        by="subcluster")
#' combined_obj <- infercnv::fusion_HMM_report_two_methods(expression_file_path = out_dir,
#'                                                         expression_file_token = "17_HMM_predHMMi6.hmm_mode-samples",
#'                                                         allele_file_path = out_dir,
#'                                                         allele_file_token = file_gene_token,
#'                                                         infercnv_obj = infercnv_object_example,
#'                                                         infercnv_allele_obj = infercnv_object_allele_gene_example,
#'                                                         method = "common",
#'                                                         type = "i6",
#'                                                         cluster_by_groups=F,
#'                                                         output_path = out_dir,
#'                                                         output_prefix = "combined")

fusion_HMM_report_two_methods <- function(expression_file_path,
                                          expression_file_token,
                                          allele_file_path,
                                          allele_file_token,
                                          # infercnv_obj,
                                          infercnv_allele_obj,
                                          method = c("union", "common"),
                                          HMM_type = c("i6", "i3"),
                                          output_path,
                                          output_prefix,
                                          HMM_report_by = 'subcluster',
                                          ...) {
  
    #allele_mode <- match.arg(allele_mode)
    method <- match.arg(method)
    HMM_type <- match.arg(HMM_type)
    
    flog.info("Loading candidate boundaries inferred from HMM ...")
    
    allele_cnv_regions_df <- read.table(file.path(allele_file_path, paste0(allele_file_token,".pred_cnv_regions.dat")),
                                        header = T,
                                        check.names = FALSE,
                                        sep="\t",
                                        stringsAsFactors = TRUE)
    expression_cnv_regions_df <- read.table(file.path(expression_file_path, paste0(expression_file_token,".pred_cnv_regions.dat")),
                                            header = T,
                                            check.names = FALSE,
                                            sep="\t",
                                            stringsAsFactors = TRUE)
    
    allele_region_gr <- GRanges(seqnames = allele_cnv_regions_df[["chr"]],
                                IRanges(as.numeric(as.character(allele_cnv_regions_df[["start"]])),
                                        as.numeric(as.character(allele_cnv_regions_df[["end"  ]]))))
    allele_region_gr$state <- allele_cnv_regions_df$state
    
    if (HMM_type == "i6") {
        allele_region_gr$state <- 2
    } else {
        allele_region_gr$state <- 1
    }
    
    allele_region_gr <- allele_region_gr %>% split(allele_cnv_regions_df$cell_group_name) %>% 
                                             GRangesList()
      
    expression_region_gr <- GRanges(seqnames = expression_cnv_regions_df[["chr"]],
                                    IRanges(as.numeric(as.character(expression_cnv_regions_df[["start"]])),
                                            as.numeric(as.character(expression_cnv_regions_df[["end"]]))))
    expression_region_gr$state <- expression_cnv_regions_df$state
    expression_region_gr <- expression_region_gr %>% split(expression_cnv_regions_df$cell_group_name) %>% 
                                                     GRangesList()
    
    union_gr <- union(names(allele_region_gr), names(expression_region_gr))
    
    if (mean(union_gr %in% names(allele_region_gr)) != 1) {
        add_gr <- setdiff(union_gr, names(allele_region_gr))
        for (i in add_gr) {
            tmp <- GRangesList(GRanges())
            names(tmp) <- i
            allele_region_gr <- c(allele_region_gr, tmp)
        }
    }
    
    if (mean(union_gr %in% names(expression_region_gr)) != 1) {
        add_gr <- setdiff(union_gr, names(expression_region_gr))
        for (i in add_gr) {
            tmp <- GRangesList(GRanges())
            names(tmp) <- i
            expression_region_gr <- c(expression_region_gr, tmp)
        }
    }
    
    flog.info("Integrating boundaries from two methods ...")
    
    combined_gr <- lapply(union_gr, function(x) {
        
        tmp_gr <- suppressWarnings(c(expression_region_gr[[x]], allele_region_gr[[x]])) %>% disjoin()
        
        tmp_gr$allele_state <- ifelse(HMM_type == "i6",
                                      3,
                                      2)
        allele_idx <- findOverlaps(tmp_gr, allele_region_gr[[x]])
        tmp_gr$allele_state[queryHits(allele_idx)] <- allele_region_gr[[x]]$state[subjectHits(allele_idx)]
        
        tmp_gr$expression_state <- ifelse(HMM_type == "i6",
                                          3,
                                          2)
        expression_idx <- findOverlaps(tmp_gr, expression_region_gr[[x]])
        tmp_gr$expression_state[queryHits(expression_idx)] <- expression_region_gr[[x]]$state[subjectHits(expression_idx)]
        
        if (HMM_type == "i6") {
            
            tmp_gr$state <- ifelse(tmp_gr$allele_state == tmp_gr$expression_state,
                                   tmp_gr$allele_state,
                                   ifelse(tmp_gr$allele_state == 2,
                                          ifelse(tmp_gr$expression_state == 3,
                                                 7,
                                                 8),
                                          tmp_gr$expression_state))
            
        } else {
            tmp_gr$state <- ifelse(tmp_gr$allele_state == tmp_gr$expression_state,
                                   tmp_gr$allele_state,
                                   ifelse(tmp_gr$allele_state == 1,
                                          ifelse(tmp_gr$expression_state == 2,
                                                 4,
                                                 5),
                                          tmp_gr$expression_state))
        }
        return(tmp_gr)
                                      
    }) %>% GRangesList()
    names(combined_gr) <- union_gr
    
    # mcmc_allele <- initialize_allele_mcmc(file_path = allele_file_path,
    #                                       file_token = allele_file_token,
    #                                       infercnv_allele_obj = infercnv_allele_obj,
    #                                       mode = allele_mode)
    
    if (method == "common") {
        infercnv_allele_obj <- remove_genes(infercnv_allele_obj, which(!rownames(infercnv_allele_obj@gene_order) %in% 
                                                                       infercnv_allele_obj@SNP_info$gene))
    }
    
    infercnv_gene_GR <- GRanges(seqnames = infercnv_allele_obj@gene_order[[C_CHR]],
                                IRanges(as.numeric(as.character(infercnv_allele_obj@gene_order[[C_START]])),
                                        as.numeric(as.character(infercnv_allele_obj@gene_order[[C_STOP ]]))))
    
    flog.info("Reassigning states estimated from combined boundaries ...")
    
    infercnv_allele_obj@expr.data[,] <- ifelse(HMM_type == "i6",
                                               3,
                                               2)
    
    for (i in names(combined_gr)) {
      
        cell_list <- infercnv_allele_obj@tumor_subclusters$subclusters %>% unlist(recursive = F)
        cell_idx <- cell_list[[i]]
        
        map_idx <- findOverlaps(infercnv_gene_GR, combined_gr[[i]])
        duplicated_idx <- which(!duplicated(queryHits(map_idx)))
        
        gene_idx <- queryHits(map_idx)[duplicated_idx]
        state_idx <- subjectHits(map_idx)[duplicated_idx]
        state_value <- combined_gr[[i]]$state[state_idx]
        
        infercnv_allele_obj@expr.data[gene_idx,cell_idx] <- matrix(state_value,
                                                                   nrow = length(state_value),
                                                                   ncol = length(cell_idx))
    }
    
    # infercnv_hmm_gene_GR <- GRanges(seqnames = infercnv_hmm_obj@gene_order[[C_CHR]],
    #                                 IRanges(as.numeric(as.character(infercnv_hmm_obj@gene_order[[C_START]])),
    #                                         as.numeric(as.character(infercnv_hmm_obj@gene_order[[C_STOP]]))))
    # 
    # for(i in seq_along(mcmc_allele@cell_gene)){
    #   #browser()
    #   cell_list <- mcmc_allele@cell_gene[[i]]$Cells
    #   gene_list <- infercnv_hmm_gene_GR %over% range(mcmc_allele@SNP_info[mcmc_allele@cell_gene[[i]]$Genes])
    #   
    #   infercnv_hmm_obj@expr.data[gene_list,cell_list] <- ifelse(type == "i6", 2, 1)
    #   
    # }
    
    infercnv::plot_cnv(infercnv_allele_obj, 
                       out_dir = output_path,
                       preset = ifelse(HMM_type == "i6",
                                         "HMM_i6_combined",
                                         "HMM_i3_combined"),
                       # x.center = ifelse(HMM_type == "i6",
                       #                   3,
                       #                   2),
                       # x.range = c(1,ifelse(HMM_type == "i6",
                       #                    8,
                       #                    5)),
                       output_format = "png",
                       output_filename = "fusion",
                       png_res = 300,
                       custom_color_pal = color.palette(c("DarkBlue", "white","DarkRed")),
                       ...)
    
    infercnv:::generate_cnv_region_reports(infercnv_allele_obj,
                                           output_filename_prefix = output_prefix,
                                           out_dir = output_path,
                                           ignore_neutral_state = ifelse(HMM_type == "i6", 3, 2),
                                           by = HMM_report_by)
    
    return(infercnv_allele_obj)
  
}

