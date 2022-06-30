#' MCMC infercnv_allele class
#' 
#' @description This class extends the functionality of infercnv_allele class
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
#' @export
#' 
MCMC_infercnv_allele <- methods::setClass(
  "MCMC_infercnv_allele", 
  slots = c(bugs_model = "character",
            cell_gene = "list",
            cnv_regions = "factor"),
  contains = "infercnv_allele")

# file_path the path to HMM reports
# file_token file name
# infercnv_allele_obj based obj
# two modes: snp and gene
initialize_allele_mcmc <- function(file_path,
                                   file_token,
                                   infercnv_allele_obj,
                                   mode = c("snp_level","gene_level")){
  
  mode <- match.arg(mode)
  
  ## Load the files for cnv predictions
  cell_groups_df <- read.table(file.path(file_path, 
                                         paste0(file_token,".cell_groupings")), 
                               header = T, check.names = FALSE, sep="\t")
  pred_cnv_genes_df <- read.table(file.path(file_path, 
                                            paste0(file_token,".pred_cnv_genes.dat")),
                                  header = T, check.names = FALSE, sep="\t", stringsAsFactors = TRUE)
  
  flog.info(sprintf("Initializing MCMC_infercnv_allele obj at %s ...", mode))
  
  mcmc_snp <- new("MCMC_infercnv_allele",
                  infercnv_allele_obj)
  mcmc_snp@bugs_model <- ifelse(mode == "snp_level",
                                system.file("BUGS_SNP_Model",package = "infercnv"),
                                system.file("BUGS_SNP2Gene_Model",package = "infercnv"))
    
  mcmc_snp@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
  
  mcmc_snp <- getGenesCells_allele(mcmc_snp,
                                   pred_cnv_genes_df, 
                                   cell_groups_df,
                                   mode = mode)
  return(mcmc_snp)
  
  # if(mode == "snp_level"){
  #   
  #   flog.info("Initializing MCMC_infercnv_allele_snp obj ...")
  #   
  #   mcmc_snp <- new("MCMC_infercnv_allele_snp",
  #                   infercnv_allele_obj)
  #   mcmc_snp@bugs_model <- system.file("BUGS_SNP_Model",package = "infercnv")
  #   mcmc_snp@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
  #   
  #   mcmc_snp <- getGenesCells_allele(mcmc_snp,
  #                                    pred_cnv_genes_df, 
  #                                    cell_groups_df,
  #                                    mode = mode)
  #   return(mcmc_snp)
  #   
  # } else if(mode == "gene_level"){
  #   flog.info("Initializing MCMC_infercnv_allele_gene obj ...")
  #   
  #   mcmc_snp_gene <- new("MCMC_infercnv_allele_gene",
  #                        infercnv_allele_obj)
  #   mcmc_snp_gene@bugs_model <- system.file("BUGS_SNP2Gene_Model",package = "infercnv")
  #   mcmc_snp_gene@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
  #   
  #   mcmc_snp_gene <- getGenesCells_allele(mcmc_snp_gene,
  #                                         pred_cnv_genes_df, 
  #                                         cell_groups_df,
  #                                         mode = mode)
  #   return(mcmc_snp_gene)
  # }
}

getGenesCells_allele <- function(obj, pred_cnv_genes_df, cell_groups_df, 
                                 mode = c("snp_level", "gene_level")){
  mode = match.arg(mode)
  
  obj@cell_gene <- lapply(obj@cnv_regions, function(x){
    #browser()
    current_cnv <- pred_cnv_genes_df[which(x == pred_cnv_genes_df$gene_region_name),]
    genes <- current_cnv$gene
    
    gene_idx <- which(row.names(obj@expr.data) %in% genes)
    sub_cells <- unique(current_cnv$cell_group_name)
    cells_idx <- which(colnames(obj@expr.data) %in% cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)
    state <- unique(current_cnv$state)
    
    if(mode == "snp_level"){
      
      ## extract affected genes
      genes.of.interest <- unique(obj@SNP_info[genes]$gene)
      
      ## associate each gene factor with a set of snps
      candidate_snp <- obj@SNP_info[genes]
      genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
        names(candidate_snp)[which(candidate_snp$gene %in% genes.of.interest[i])]
      })
      names(genes2snps.dict) <- genes.of.interest
      
      multi_arrays <- set_array(r.maf = obj@count.data[gene_idx,cells_idx], 
                                n.sc = obj@coverage.data[gene_idx,cells_idx],
                                genes2snps.dict = genes2snps.dict, 
                                numCells = length(cells_idx))
      
      return(list("cnv_regions" = x, 
                  "SNPs" = gene_idx, 
                  "Cells" = cells_idx, 
                  "State" = state,
                  
                  "r.array" = multi_arrays$r.array,
                  "n.sc.array" = multi_arrays$n.sc.array,
                  "l.array" = multi_arrays$l.array,
                  "n.bulk.array" = multi_arrays$n.bulk.array,
                  "I.j"=multi_arrays$I.j))
      
    } else{
      return(list("cnv_regions" = x, 
                  "Genes" = gene_idx, 
                  "Cells" = cells_idx, 
                  "State" = state))
    }
  })
  return(obj)
}

## the idea of array initialization is borrowed from HoneyBadger
set_array <- function(r.maf, n.sc, genes2snps.dict, numCells){
  #browser()
  I.j <- unlist(lapply(genes2snps.dict, length))
  numGenes <- length(genes2snps.dict)
  numSnpsPerGene <- max(I.j)
  l.maf <- rowSums(r.maf > 0)
  n.bulk <- rowSums(n.sc > 0)

  ## j, i, k
  r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      r.array[i,s,] <- r.maf[snpst[s],]
    }
  }
  n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      n.sc.array[i,s,] <- n.sc[snpst[s],]
    }
  }
  l.array <- array(0, c(numGenes, numSnpsPerGene))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      l.array[i,s] <- l.maf[snpst[s]]
    }
  }
  n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      n.bulk.array[i,s] <- n.bulk[snpst[s]]
    }
  }
  return(list("r.array"=r.array,
              "n.sc.array"=n.sc.array,
              "l.array"=l.array,
              "n.bulk.array"=n.bulk.array,
              "I.j"=I.j))
}

#' @title inferCNVAlleleBayesNet: Run Bayesian Network Mixture Model Leveraging Allele Data 
#' To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of two Copy Number Variation states 
#' (i2 states: 1,Deletion; 2,Neutral) for CNV's identified by inferCNV's HMM. 
#' Posterior probabilities are found for the entire CNV cluster and each individual cell line in the CNV.
#' 
#' @param file_path Location of the directory of the inferCNV_allele outputs.
#' 
#' @param file_token (string) String token that contains some info on settings used to name allele-based files.
#' 
#' @param infercnv_allele_obj InferCNV allele object.
#' 
#' @param allele_mode The type of allele data provided. "snp_level" or "gene_level".
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param cores Option to run parallel by specifying the number of cores to be used. (Default: 5)
#'
#' @export

inferCNVAlleleBayesNet <- function(file_path,
                                   file_token,
                                   infercnv_allele_obj,
                                   allele_mode = c("snp_level","gene_level"),
                                   output_path,
                                   cores = 5){
  
  if(!dir.exists(file.path(output_path))){
    dir.create(file.path(output_path), recursive = T)
    flog.info(paste("Creating the following Directory:", output_path))
  }
  
  allele_mode <- match.arg(allele_mode)
  
  mcmc_allele <- initialize_allele_mcmc(file_path = file_path,
                                        file_token = file_token,
                                        infercnv_allele_obj = infercnv_allele_obj,
                                        mode = allele_mode)
  flog.info(paste("The number of affected regions:", length(mcmc_allele@cell_gene)))
  
  flog.info("Start running Gibbs sampling ")
  
  start_time <- Sys.time()
  if(allele_mode == "snp_level"){
    
    mclapply(mcmc_allele@cell_gene, function(x){
      #browser()
      samples <- run_allele_snp_mcmc(bugs = mcmc_allele@bugs_model,
                                     r.array = x$r.array,
                                     n.sc.array = x$n.sc.array,
                                     l.array = x$l.array,
                                     n.bulk.array = x$n.bulk.array,
                                     nGenes = length(x$I.j),
                                     nCells = length(x$Cells),
                                     I.j = x$I.j,
                                     pe = 0.1,
                                     mono = 0.7)
      plot_mcmc(samples = samples,
                region = x$cnv_regions,
                output_path = output_path)
    }, mc.cores = cores)
  } else{
    mclapply(mcmc_allele@cell_gene, function(x){
      #browser()
      samples <- run_allele_gene_mcmc(bugs = mcmc_allele@bugs_model,
                                      r = mcmc_allele@count.data[x$Genes,
                                                                 x$Cells,
                                                                 drop = F],
                                      n.sc = mcmc_allele@coverage.data[x$Genes,
                                                                       x$Cells,
                                                                       drop = F],
                                      nGenes = length(x$Genes),
                                      nCells = length(x$Cells),
                                      pe = 0.1,
                                      mono = 0.7)
      plot_mcmc(samples = samples,
                region = x$cnv_regions,
                output_path = output_path)
    }, mc.cores = cores)
  }
  end_time <- Sys.time()
  flog.info(paste0("MCMC running time: ",
                   difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
}

run_allele_snp_mcmc <- function(bugs,
                                r.array,
                                n.sc.array,
                                l.array,
                                n.bulk.array,
                                nGenes,
                                nCells,
                                I.j,
                                pe = 0.1,
                                mono = 0.7){
  #browser()
  
  input <- list(
    'r' = r.array,
    'l' = l.array,
    'n.sc' = n.sc.array,
    'n.bulk' = n.bulk.array,
  
    'J' = nGenes,  # how many genes
    'K' = nCells,  # how many cells
    'I.j' = I.j,
    'pseudo' = pe,
    'mono' = mono)
  
  model <- rjags::jags.model(bugs,
                             data = input,
                             n.chains = 3, n.adapt = 500,
                             quiet = F)
  update(model, 200, quiet = F)
  parameters <- c('theta', 'epsilon')
  samples <- rjags::coda.samples(model, parameters, 
                                 n.iter=1e3, 
                                 thin = 10,
                                 quiet = F)
  
  return(samples)
}

run_allele_gene_mcmc <- function(bugs,
                                 r,
                                 n.sc,
                                 nGenes,
                                 nCells,
                                 pe = 0.1,
                                 mono = 0.7){
  #browser()
  input <- list(
    'r' = r,
    'n.sc' = n.sc,
    'l' = rowSums(r > 0),
    'n.bulk' = rowSums(n.sc > 0),
    'J' = nGenes,  # how many genes (snp)
    'K' = nCells, # how many cells
    'pseudo' = pe,
    'mono' = mono)
  
  model <- rjags::jags.model(bugs,
                             data = input,
                             n.chains = 3, n.adapt = 500,
                             quiet = F)
  update(model, 200, quiet = F)
  parameters <- c('theta', 'epsilon')
  samples <- rjags::coda.samples(model, parameters, 
                                 n.iter=1e3, 
                                 thin = 10,
                                 quiet = F)
  
  return(samples)
  
}

plot_mcmc <- function(samples,
                      region,
                      output_path){
  #browser()
  theta_samples <- lapply(samples, function(x){
    tmp <- x[,grepl("theta", colnames(x))]
  })
  
  theta_samples_df <- do.call(rbind, theta_samples) #%>% reshape::melt()
  colnames(theta_samples_df) <- seq_len(ncol(theta_samples_df))
  # if(ncol(theta_samples_df) == 3){
  #   colnames(theta_samples_df) <- c("del", "neutral","amp")
  # } else if(ncol(theta_samples_df) == 7){
  #   colnames(theta_samples_df) <- c("del x2", "del x1",
  #                                   "neutral",
  #                                   "amp x1", "amp x2", "amp x3",
  #                                   "cnLOH")
  # } else if(ncol(theta_samples_df) == 4){
  #   colnames(theta_samples_df) <- c("del", 
  #                                   "neutral",
  #                                   "amp", 
  #                                   "cnLOH")
  #   
  # } else if(ncol(theta_samples_df) == 6){
  #   colnames(theta_samples_df) <- c("del x2", "del x1",
  #                                   "neutral",
  #                                   "amp x1", "amp x2", "amp x3")
  # } else{colnames(theta_samples_df) <- c("del", "neutral")}
  #    
  theta_samples_df <- theta_samples_df %>% reshape::melt()
  colnames(theta_samples_df) <- c("index", "State", "Probability")
  # if(length(unique(theta_samples_df$State)) == 3){
  #   theta_samples_df$State <- factor(theta_samples_df$State,
  #                                    levels = c("del", "neutral","amp"))
  # } else if(length(unique(theta_samples_df$State)) == 7){
  #   theta_samples_df$State <- factor(theta_samples_df$State,
  #                                    levels = c("del x2", "del x1", 
  #                                               "neutral",
  #                                               "amp x1", "amp x2", "amp x3",
  #                                               "cnLOH"))
  # } else if(length(unique(theta_samples_df$State)) == 4){
  #   theta_samples_df$State <- factor(theta_samples_df$State,
  #                                    levels = c("del", "neutral","amp","cnLOH"))
  # } else if(length(unique(theta_samples_df$State)) == 6){
  #   theta_samples_df$State <- factor(theta_samples_df$State,
  #                                    levels = c("del x2", "del x1", 
  #                                               "neutral",
  #                                               "amp x1", "amp x2", "amp x3"))
  # } else{
  #   theta_samples_df$State <- factor(theta_samples_df$State,
  #                                    levels = c("del", "neutral"))
  # }
  # 
  # stat_box_data <- function(y, upper_limit = max(theta_samples_df$Prob) * 1.15) {
  #   return( 
  #     data.frame(
  #       y = 0.95 * upper_limit,
  #       label = paste('count =', length(y), '\n',
  #                     'mean =', round(mean(y), 3), '\n')
  #     )
  #   )
  # }
  
  pdf(file = file.path(output_path,
                       paste0(region,"_DiagnosticPlots.pdf")), onefile = TRUE)
  coda::autocorr.plot(coda::as.mcmc.list(theta_samples))
  coda::gelman.plot(coda::mcmc.list(theta_samples))
  dev.off()
  
  pdf(file = file.path(output_path,
                       paste0(region,"_DensityPlots.pdf")), onefile = TRUE)
  plot(samples)
  dev.off()
  
  # prob_plot <- ggplot(theta_samples_df,
  #                     aes(y=Prob, x=State, color = State)) +
  #   geom_boxplot() +
  #   #ylim(c(0,1)) +
  #   stat_summary(fun.data = stat_box_data, 
  #                geom = "text", 
  #                hjust = 0.5,
  #                vjust = 0.9) +
  #   theme_bw() + 
  #   ggtitle(region)
  states <- as.factor(theta_samples_df$State)
  prob_plot <- ggplot(data = theta_samples_df, aes_string(y = 'Probability', x= 'State', fill = 'states')) +
    geom_boxplot()+
    labs(title = region) +
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(file = file.path(output_path,
                       paste0(region,"_CNVProbPlots.pdf")), onefile = TRUE)
  print(prob_plot)
  dev.off()
  # ggsave(file.path(output_path,
  #                  paste0(region,"_CNVProbPlots.pdf")),
  #        prob_plot)
}

