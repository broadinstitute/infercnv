#' MCMC infercnv_combined class
#' 
#' @description This class extends the functionality of infercnv/infercnv_allele class
#' aiming to cooperate MCMC object leveraging Bayesian model
#' 
#' Slots in the MCMC_infercnv_allele object (besides infercnv) include:
#' 
#' @slot bugs_model a string path to where the bug file locates.
#' 
#' @slot cell_gene cell_gene.
#' 
#' @slot cnv_regions cnv_regions.
#' 
#' @slot infercnv_obj gene_obj.
#' 
#' @slot infercnv_obj_mu infercnv_obj_mu.
#' 
#' @slot infercnv_obj_sig infercnv_obj_sig.
#' 
#' @slot infercnv_allele_obj allele_obj.
#' 
#' @export
#' 
MCMC_infercnv_combined <- methods::setClass(
  "MCMC_infercnv_combined", 
  slots = c(bugs_model = "character",
            cell_gene = "list",
            cnv_regions = "factor",
            infercnv_obj = "infercnv",
            infercnv_obj_mu = "numeric",
            infercnv_obj_sig = "numeric",
            infercnv_allele_obj = "infercnv_allele"))

# file_path the path to HMM reports
# file_token file name
# infercnv based obj
# infercnv_allele based obj
# two modes: snp and gene
initialize_combined_mcmc <- function(file_path,
                                     file_token,
                                     infercnv_obj,
                                     infercnv_obj_mu,
                                     infercnv_obj_sig,
                                     infercnv_allele_obj,
                                     mode = c("snp_level","gene_level"),
                                     type = c("i6", "i3"),
                                     enable_cnLOH = TRUE){
  
  mode <- match.arg(mode)
  type <- match.arg(type)
  
  ## Load the files for cnv predictions
  cell_groups_df <- read.table(file.path(file_path, 
                                         paste0(file_token,".cell_groupings")), 
                               header = T, check.names = FALSE, sep="\t")
  pred_cnv_genes_df <- read.table(file.path(file_path, 
                                            paste0(file_token,".pred_cnv_genes.dat")),
                                  header = T, check.names = FALSE, sep="\t", stringsAsFactors = TRUE)
  
  flog.info(sprintf("Initializing MCMC_infercnv_combined obj at %s ...", mode))
  
  mcmc_combined <- new("MCMC_infercnv_combined",
                        infercnv_obj = infercnv_obj,
                        infercnv_obj_mu = infercnv_obj_mu,
                        infercnv_obj_sig = infercnv_obj_sig,
                        infercnv_allele_obj = infercnv_allele_obj)
  
  mcmc_combined@bugs_model <- ifelse(type == "i6",
                                     ifelse(mode == "snp_level",
                                            ifelse(enable_cnLOH,
                                                   system.file("BUGS_Combined_Model_i7",package = "infercnv"),
                                                   system.file("BUGS_Combined_Model_i6",package = "infercnv")),
                                            ifelse(enable_cnLOH,
                                                   system.file("BUGS_SNP2Gene_Combined_Model_i7",package = "infercnv"),
                                                   system.file("BUGS_SNP2Gene_Combined_Model_i6",package = "infercnv"))),
                                     ifelse(mode == "snp_level",
                                            ifelse(enable_cnLOH,
                                                   system.file("BUGS_Combined_Model_i4",package = "infercnv"),
                                                   system.file("BUGS_Combined_Model_i3",package = "infercnv")),
                                            ifelse(enable_cnLOH,
                                                   system.file("BUGS_SNP2Gene_Combined_Model_i4",package = "infercnv"),
                                                   system.file("BUGS_SNP2Gene_Combined_Model_i3",package = "infercnv"))))

  mcmc_combined@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
  
  mcmc_combined <- getGenesCells_combined(mcmc_combined,
                                          pred_cnv_genes_df, 
                                          cell_groups_df,
                                          mode = mode)
                                              
  return(mcmc_combined)
  
  # if(mode == "snp_level"){
  #   
  #   flog.info("initializing MCMC_infercnv_combined obj at snp level ...")
  #   
  #   mcmc_combined_snp <- new("MCMC_infercnv_combined",
  #                            infercnv_obj = infercnv_obj,
  #                            infercnv_obj_mu = infercnv_obj_mu,
  #                            infercnv_obj_sig = infercnv_obj_sig,
  #                            infercnv_allele_obj = infercnv_allele_obj)
  #   mcmc_combined_snp@bugs_model <- ifelse(type == "i6",
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
  #   mcmc_combined_gene@bugs_model <- ifelse(type == "i6",
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

getGenesCells_combined <- function(obj, pred_cnv_genes_df, cell_groups_df, 
                                   mode = c("snp_level", "gene_level")){
  mode = match.arg(mode)
  
  obj@cell_gene <- lapply(obj@cnv_regions, function(x){
    #browser()
    current_cnv <- pred_cnv_genes_df[which(x == pred_cnv_genes_df$gene_region_name),]
    genes <- current_cnv$gene
    
    infercnv_gene_idx <- which(row.names(obj@infercnv_obj@expr.data) %in% genes)
    
    if(mode == "snp_level"){
      infercnv_allele_gene <- names(obj@infercnv_allele_obj@SNP_info)[obj@infercnv_allele_obj@SNP_info$gene %in% genes]
      infercnv_allele_gene_idx <- which(obj@infercnv_allele_obj@SNP_info$gene %in% 
                                          infercnv_allele_gene)
    } else{
      infercnv_allele_gene_idx <- which(obj@infercnv_allele_obj@SNP_info$gene %in% genes)
    }
    
    
    sub_cells <- unique(current_cnv$cell_group_name)
    infercnv_cells_idx <- which(colnames(obj@infercnv_obj@expr.data) %in% 
                                  cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
    infercnv_allele_cells_idx <- which(colnames(obj@infercnv_allele_obj@expr.data) %in% 
                                         cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
    state <- unique(current_cnv$state)
    
    if(length(infercnv_allele_gene_idx) == 0){
      
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
    
    if(mode == "snp_level"){
      
      ## extract affected genes
      genes.of.interest <- unique(obj@infercnv_allele_obj@SNP_info[infercnv_allele_gene_idx]$gene)
      
      ## associate each gene factor with a set of snps
      candidate_snp <- obj@infercnv_allele_obj@SNP_info[infercnv_allele_gene_idx]
      genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
        names(candidate_snp)[which(candidate_snp$gene %in% genes.of.interest[i])]
      })
      names(genes2snps.dict) <- genes.of.interest
      
      multi_arrays <- set_array(r.maf = obj@infercnv_allele_obj@count.data[infercnv_allele_gene_idx,
                                                                           infercnv_allele_cells_idx,
                                                                           drop = F], 
                                n.sc = obj@infercnv_allele_obj@coverage.data[infercnv_allele_gene_idx,
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
      
    } else{
      return(list("cnv_regions" = x, 
                  
                  "infercnv_Genes" = infercnv_gene_idx,
                  "infercnv_Cells" = infercnv_cells_idx,
                  
                  "infercnv_allele_Genes" = infercnv_allele_gene_idx, 
                  "infercnv_allele_Cells" = infercnv_allele_cells_idx, 
                  
                  "State" = state))
    }
  })
  return(obj)
}

#' @title inferCNVCombinedBayesNet: Run Bayesian Network Mixture Model Leveraging both
#' Expression and Allele Data To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of two Copy Number Variation states 
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
#' @param type The type of combined model running in Bayesian Method. i6 or i3.
#' 
#' @param enable_cnLOH Option for detecting cnLOH events by adding one more state in Bayesian Method. (Default: True)
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param cores Option to run parallel by specifying the number of cores to be used. (Default: 5)
#'
#' @export

inferCNVCombinedBayesNet <- function(combined_file_path,
                                     combined_file_token,
                                     infercnv_obj,
                                     infercnv_allele_obj,
                                     allele_mode = c("snp_level","gene_level"),
                                     type = c("i6", "i3"),
                                     enable_cnLOH = TRUE,
                                     output_path,
                                     cores = 5){
  
  allele_mode <- match.arg(allele_mode)
  type <- match.arg(type)
  
  if(!dir.exists(file.path(output_path))){
    dir.create(file.path(output_path), recursive = T)
    flog.info(paste("Creating the following Directory:", output_path))
  }
  
  if(type == "i6"){
    cnv_mean_sd = infercnv:::get_spike_dists(infercnv_obj@.hspike)
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
  } else{
    cnv_mean_sd = infercnv:::.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj)
    mean <- c(cnv_mean_sd$mu - cnv_mean_sd$mean_delta,
              cnv_mean_sd$mu,
              cnv_mean_sd$mu + cnv_mean_sd$mean_delta)
    
    sd <- c(1/(cnv_mean_sd$sigma^2),
            1/(cnv_mean_sd$sigma^2),
            1/(cnv_mean_sd$sigma^2))
  }
  
  mcmc_combined <- initialize_combined_mcmc(file_path = combined_file_path,
                                            file_token = combined_file_token,
                                            infercnv_obj = infercnv_obj,
                                            infercnv_obj_mu = mean,
                                            infercnv_obj_sig = sd,
                                            infercnv_allele_obj = infercnv_allele_obj,
                                            mode = allele_mode,
                                            type = type,
                                            enable_cnLOH = enable_cnLOH)
                                                       
  
  flog.info(paste("The number of affected regions:", length(mcmc_combined@cell_gene)))
  
  flog.info("Start running Gibbs sampling ")
  
  start_time <- Sys.time()
  
  mclapply(mcmc_combined@cell_gene, function(x){
    #browser()
    if(allele_mode == "snp_level"){
      # if(is.null(x$r.array)){
      #   mcmc_combined@bugs_model <- ifelse(type == "i6",
      #                                      system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
      #                                      system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv"))
      # } 
      samples <- run_combined_snp_mcmc(bugs = ifelse(is.null(x$r.array),
                                                     ifelse(type == "i6",
                                                            system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
                                                            system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv")),
                                                     mcmc_combined@bugs_model),
                                       type = type,
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
                                       gexp = mcmc_combined@infercnv_obj@expr.data[x$infercnv_Genes,
                                                                                   x$infercnv_Cells,
                                                                                   drop = F],
                                       mu = mcmc_combined@infercnv_obj_mu,
                                       sig = mcmc_combined@infercnv_obj_sig,
                                       pe = 0.1,
                                       mono = 0.7)
      
      plot_mcmc(samples = samples,
                region = x$cnv_regions,
                output_path = output_path)
      
    } else{
      
      # if(is.null(x$infercnv_allele_Genes)){
      #   mcmc_combined@bugs_model <- ifelse(type == "i6",
      #                                      system.file("BUGS_Mixture_Model_i6_gene_test",package = "infercnv"),
      #                                      system.file("BUGS_Mixture_Model_i3_gene_test",package = "infercnv"))
      # } 
      samples <- run_combined_gene_mcmc(bugs = ifelse(is.null(x$infercnv_allele_Genes),
                                                      ifelse(type == "i6",
                                                             system.file("BUGS_Mixture_Model_i6_gene",package = "infercnv"),
                                                             system.file("BUGS_Mixture_Model_i3_gene",package = "infercnv")),
                                                      mcmc_combined@bugs_model),
                                        type = type,
                                        enable_cnLOH = ifelse(is.null(x$infercnv_allele_Genes),
                                                              FALSE,
                                                              enable_cnLOH),
                                        r = mcmc_combined@infercnv_allele_obj@count.data[x$infercnv_allele_Genes,
                                                                                         x$infercnv_allele_Cells,
                                                                                         drop = F],
                                        n.sc = mcmc_combined@infercnv_allele_obj@coverage.data[x$infercnv_allele_Genes,
                                                                                               x$infercnv_allele_Cells,
                                                                                               drop = F],
                                        nAlleleGenes = length(x$infercnv_allele_Genes),
                                        nGenes = length(x$infercnv_Genes),
                                        nCells = length(x$infercnv_Cells),
                                        gexp = mcmc_combined@infercnv_obj@expr.data[x$infercnv_Genes,
                                                                                    x$infercnv_Cells,
                                                                                    drop = F],
                                        mu = mcmc_combined@infercnv_obj_mu,
                                        sig = mcmc_combined@infercnv_obj_sig,
                                        pe = 0.1,
                                        mono = 0.7)
      
      plot_mcmc(samples = samples,
                region = x$cnv_regions,
                output_path = output_path)
    }
  }, mc.cores = cores)
  #})
  end_time <- Sys.time()
  flog.info(paste0("MCMC running time: ",
                   difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
}

run_combined_snp_mcmc <- function(bugs,
                                  type = c("i6","i3"),
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
                                  mono = 0.7){
  
  type <- match.arg(type)
  
  #browser()
  input <- list(
    'r' = r.array,
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
  
  if(type == "i6"){
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
  samples <- rjags::coda.samples(model, parameters, 
                                 n.iter=1e3, 
                                 thin = 10,
                                 quiet = F)
  
  return(samples)
}

run_combined_gene_mcmc <- function(bugs,
                                   type = c("i6","i3"),
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
                                   mono = 0.7){
  
  type <- match.arg(type)
  
  #browser()
  input <- list(
    'r' = r,
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
  
  if(type == "i6"){
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
  
  samples <- rjags::coda.samples(model, parameters, 
                                 n.iter=1e3, 
                                 thin = 10,
                                 quiet = F)
  
  return(samples)
  
}

#' @title fusion_HMM_report_two_methods 
#'
#' @description 
#' 
#' @param infercnv_hmm_obj InferCNV HMM object.
#' 
#' @param allele_file_path Location of the directory of the inferCNV_allele outputs.
#' 
#' @param allele_file_token (string) String token that contains some info on settings used to name allele-based files.
#' 
#' @param infercnv_allele_obj InferCNV allele object.
#' 
#' @param allele_mode The type of allele data provided. "snp_level" or "gene_level".
#' 
#' @param method The way integrating HMM boundaries from two methods: union/common.
#' "union" keeps those expression-specific gene that are not found in allele data during modeling,
#' while "common" removes those ones.
#' 
#' @param type The type of combined model running in Bayesian Method. i6 or i3.
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param output_prefix name prefix when naming output files.
#' 
#' @param HMM_report_by cell, consensus, subcluster (default: subcluster)
#'
#' @export

fusion_HMM_report_two_methods <- function(infercnv_hmm_obj,
                                          allele_file_path,
                                          allele_file_token,
                                          infercnv_allele_obj,
                                          allele_mode = c("snp_level","gene_level"),
                                          method = c("union", "common"),
                                          type = c("i6", "i3"),
                                          output_path,
                                          output_prefix,
                                          HMM_report_by = 'subcluster',
                                          ...){
  
  allele_mode <- match.arg(allele_mode)
  method <- match.arg(method)
  type <- match.arg(type)
  
  mcmc_allele <- initialize_allele_mcmc(file_path = allele_file_path,
                                        file_token = allele_file_token,
                                        infercnv_allele_obj = infercnv_allele_obj,
                                        mode = allele_mode)
  
  if(method == "common"){
    infercnv_hmm_obj <- remove_genes(infercnv_hmm_obj,
                                      which(! rownames(infercnv_hmm_obj@gene_order) %in% 
                                              mcmc_allele@SNP_info$gene))
  }
  
  infercnv_hmm_gene_GR <- GRanges(seqnames = infercnv_hmm_obj@gene_order[[C_CHR]],
                                  IRanges(as.numeric(as.character(infercnv_hmm_obj@gene_order[[C_START]])),
                                          as.numeric(as.character(infercnv_hmm_obj@gene_order[[C_STOP]]))))
  
  for(i in seq_along(mcmc_allele@cell_gene)){
    #browser()
    cell_list <- mcmc_allele@cell_gene[[i]]$Cells
    gene_list <- infercnv_hmm_gene_GR %over% range(mcmc_allele@SNP_info[mcmc_allele@cell_gene[[i]]$Genes])
    
    infercnv_hmm_obj@expr.data[gene_list,cell_list] <- ifelse(type == "i6", 2, 1)
    
  }
  
  infercnv::plot_cnv(infercnv_hmm_obj, 
                     out_dir=output_path,
                     ...)
  
  infercnv:::generate_cnv_region_reports(infercnv_hmm_obj,
                                         output_filename_prefix=output_prefix,
                                         out_dir=output_path,
                                         ignore_neutral_state = ifelse(type == "i6", 3, 2),
                                         by=HMM_report_by)
  
  return(infercnv_hmm_obj)
  
}
