#' MCMC infercnv_allele class
#' 
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of two Copy Number Variation states 
#' (i2 states: 1,Deletion; 2,Neutral) for CNV's identified by inferCNV_allele's HMM.
#' 
#' @slot bugs_model BUGS model.
#' 
#' @slot cell_gene List containing the Cells and Genes that make up each CNV.
#' 
#' @slot cnv_probabilities Probabilities of each CNV belonging to a particular state from 0 (least likely) to 1 (most likely).
#' 
#' @slot cell_probabilities Probabilities of each cell being in a particular state, from 0 (least likely) to 1 (most likely).
#' 
#' @slot cnv_regions ID for each CNV found by the HMM.
#' 
#' @export
#' 
MCMC_infercnv_allele <- methods::setClass(
  "MCMC_infercnv_allele", 
  slots = c(bugs_model = "character",
            cell_gene = "list",
            cnv_probabilities = "list",
            cell_probabilities = "list",
            cnv_regions = "factor"),
  contains = "infercnv_allele")

# file_path the path to cnv_reports
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

# Function for each individule cell probabilities, marginalize over the EPSILONS -- my version
my_cell_prob <- function(combined_samples) {
  epsilons <- combined_samples[,grepl('epsilon', colnames(combined_samples))]
  state <- combined_samples[,grepl('theta', colnames(combined_samples))]
  #print(paste("Epsilons: ", dim(epsilons)))
  epsilon_state_frequencies <- apply(as.data.frame(epsilons), 2, function(x) table(factor(x, levels = seq_len(ncol(state)))))
  cell_probs <- epsilon_state_frequencies/colSums(epsilon_state_frequencies)
  return(cell_probs)
}

get_posterior_prob <- function(obj, mcmc){
  
  cnv_probabilities <- list()
  ## List for combining the chains in each simulation
  combined_mcmc <- list()
  ## list holding the frequency of epsilon values for each cell line
  ##  for each cnv region and subgroup
  cell_probabilities <- list()
  
  combinedMCMC <-
    for(j in seq_along(mcmc)){
      # combine the chains
      combined_mcmc[[j]] <- do.call(rbind, mcmc[[j]])
      # run function to get probabilities
      ## Thetas
      cnv_probabilities[[j]] <- cnv_prob(combined_mcmc[[j]])
      ## Epsilons
      cell_probabilities[[j]] <- my_cell_prob(combined_mcmc[[j]])
    }
  obj@cnv_probabilities <- cnv_probabilities
  obj@cell_probabilities <- cell_probabilities
  
  return(obj)
}

plot_posterior_prob <- function(obj, output_path){
  
  pdf(file = file.path(output_path,"cnvProbs.pdf"), onefile = TRUE)
  lapply(seq_along(obj@cnv_probabilities), function(i){
    print(my_plot_cnv_prob(obj@cnv_probabilities[[i]], as.character(obj@cell_gene[[i]]$cnv_regions)))
  })
  dev.off()
  
  pdf(file = file.path(output_path,"cellProbs.pdf"), onefile = TRUE)
  lapply(seq_along(obj@cell_probabilities), function(i){
    print(my_plot_cell_prob(as.data.frame(obj@cell_probabilities[[i]]), as.character(obj@cell_gene[[i]]$cnv_regions)))
  })
  dev.off()
  
}

## Fucntion to Plot the probability of each state for a CNV -- my version
my_plot_cnv_prob <- function(df, title){
  colnames(df) <- seq_len(ncol(df))
  df <- melt(df)
  colnames(df) <- c("row", "State", "Probability")
  states <- as.factor(df$State)
  ggplot2::ggplot(data = df, ggplot2::aes_string(y = 'Probability', x= 'State', fill = 'states')) +
    ggplot2::geom_boxplot()+
    ggplot2::labs(title = title) +
    ggplot2::theme(plot.title = element_text(hjust = 0.5))
}

# Function to plot the probability for each cell line of being in a particular state -- my version
my_plot_cell_prob <- function(df, title){
  df$mag = seq_len(nrow(df))
  long_data <- reshape::melt(df, id = "mag")
  long_data$mag <- as.factor(long_data$mag)
  ggplot2::ggplot(long_data, ggplot2::aes_string(x = 'variable', y = 'value', fill = 'mag'))+
    ggplot2::geom_bar(stat="identity", width = 1) +
    ggplot2::coord_flip() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),panel.border = ggplot2::element_blank(),
      axis.text=ggplot2::element_text(size=20),
      plot.title = ggplot2::element_text(hjust = 0.5,size = 22),
      #legend.position = "none",
      legend.position="bottom",
      axis.text.x = ggplot2::element_text(size = 16),
      axis.text.y = ggplot2::element_text(size = 16),
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18))+
    ggplot2::labs(title = title) +
    #fill = "CNV States") +
    ggplot2::xlab("Cell") +
    ggplot2::ylab("Probability")+
    ggplot2::labs(fill = "States")+
    ggplot2::scale_x_discrete(breaks =seq(1, ncol(df), 9))
}

plot_Diagnostics <- function(mcmc, output_path){
  
  cnvMCMCList <- lapply(seq_along(mcmc), function(i){
    lapply(mcmc[[i]], cnv_prob)
  })
  pdf(file = file.path(output_path,"CNVDiagnosticPlots.pdf"), onefile = TRUE)
  lapply(seq_along(cnvMCMCList), function(i){
    plot(coda::mcmc.list(cnvMCMCList[[i]]))
  })
  dev.off()
  
  cellProb <- function(samples) {
    epsilons <- samples[,grepl('epsilon', colnames(samples))]
  }
  cellMCMCList <- lapply(seq_along(mcmc), function(i){
    lapply(mcmc[[i]], cellProb)
  })
  pdf(file = file.path(output_path,"CellDiagnosticPlots.pdf"), onefile = TRUE)
  lapply(seq_along(cellMCMCList), function(i){
    plot(coda::mcmc.list(cellMCMCList[[i]]))
  })
  dev.off()
  
  pdf(file = file.path(output_path,"CNVautocorrelationPlots.pdf"), onefile = TRUE)
  lapply(seq_along(cnvMCMCList), function(i){
    coda::autocorr.plot(coda::mcmc.list(cnvMCMCList[[i]]))
  })
  dev.off()
  
  pdf(file = file.path(output_path,"CNVGelmanPlots.pdf"), onefile = TRUE)
  lapply(seq_along(cnvMCMCList), function(i){
    coda::gelman.plot(coda::mcmc.list(cnvMCMCList[[i]]))
  })
  dev.off()
  
}

#' @title inferCNVAlleleBayesNet: Run Bayesian Network Mixture Model Leveraging Allele Data 
#' To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of two Copy Number Variation states 
#' (i2 states: 1,Deletion; 2,Neutral) for CNV's identified by inferCNV_allele's HMM. 
#' Posterior probabilities are found for the entire CNV cluster and each individual cell line in the CNV.
#' 
#' @param file_path Location of the directory of the inferCNV_allele outputs.
#' 
#' @param file_token (string) String token that contains some info on settings used to name allele-based files.
#' 
#' @param infercnv_allele_obj InferCNV_allele object.
#' 
#' @param allele_mode The type of allele data provided. "snp_level" or "gene_level".
#' 
#' @param output_path (string) Path to where the output file should be saved to.
#' 
#' @param cores Option to run parallel by specifying the number of cores to be used. (Default: 5)
#' 
#' @return Returns a MCMC_inferCNV_allele_obj and posterior probability of being in one of two Copy Number Variation states
#' (i2 states: 1,Deletion; 2,Neutral) for CNV's identified by inferCNV_allele's HMM.
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
    
    mcmc <- mclapply(mcmc_allele@cell_gene, function(x){
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
      # plot_mcmc(samples = samples,
      #           region = x$cnv_regions,
      #           output_path = output_path)
    }, mc.cores = cores)
  } else{
    mcmc <- mclapply(mcmc_allele@cell_gene, function(x){
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
      # plot_mcmc(samples = samples,
      #           region = x$cnv_regions,
      #           output_path = output_path)
    }, mc.cores = cores)
  }
  end_time <- Sys.time()
  flog.info(paste0("MCMC running time: ",
                   difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
  
  flog.info("Retrieving the posterior probabilities of CNVs ")
  mcmc_allele <- get_posterior_prob(obj = mcmc_allele, mcmc = mcmc)
  
  saveRDS(mcmc_allele, file = file.path(output_path, "mcmc_allele.rds"))
  
  flog.info("Plotting the distribution of posterior probabilities")
  plot_posterior_prob(obj = mcmc_allele, output_path)
  
  flog.info("Plotting diagnostic statistics")
  plot_Diagnostics(mcmc = mcmc, output_path)
  
  return(mcmc_allele)
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

