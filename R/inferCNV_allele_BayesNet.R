#' MCMC infercnv_allele_snp class -- under development
#' 
#' @description This class extends the functionality of infercnv_allele class
#' aiming to coorporate MCMC object leveraging Bayesian model
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
MCMC_infercnv_allele_snp <- methods::setClass(
  "MCMC_infercnv_allele_snp", 
  slots = c(bugs_model = "character",
            cell_gene = "list",
            cnv_regions = "factor"),
  contains = "infercnv_allele")

#' MCMC infercnv_allele_gene class -- under development
#' 
#' @description This class extends the functionality of infercnv_allele class
#' aiming to coorporate MCMC object leveraging Bayesian model
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
MCMC_infercnv_allele_gene <- methods::setClass(
  "MCMC_infercnv_allele_gene", 
  slots = c(bugs_model = "character",
            cell_gene = "list",
            cnv_regions = "factor"),
  contains = "infercnv_allele")

# file_path the path to HMM reports
# file_token file name
# allele_mcmc_allele based obj
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
  
  if(mode == "snp_level"){
    
    flog.info("initializing MCMC_infercnv_allele_snp obj ...")
    
    mcmc_snp <- new("MCMC_infercnv_allele_snp",
                    infercnv_allele_obj)
    mcmc_snp@bugs_model <- system.file("BUGS_SNP_Model",package = "infercnv")
    mcmc_snp@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
    
    mcmc_snp <- getGenesCells_allele(mcmc_snp,
                                     pred_cnv_genes_df, 
                                     cell_groups_df,
                                     mode = mode)
    return(mcmc_snp)
    
  } else if(mode == "gene_level"){
    flog.info("initializing MCMC_infercnv_allele_gene obj ...")
    
    mcmc_snp_gene <- new("MCMC_infercnv_allele_gene",
                         infercnv_allele_obj)
    mcmc_snp_gene@bugs_model <- system.file("BUGS_SNP2Gene_Model",package = "infercnv")
    mcmc_snp_gene@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
    
    mcmc_snp_gene <- getGenesCells_allele(mcmc_snp_gene,
                                          pred_cnv_genes_df, 
                                          cell_groups_df,
                                          mode = mode)
    return(mcmc_snp_gene)
  }
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

## mcmc_allele mcmc_allele based obj
## output_path output path to where you get your outputs
## cores the number of cores used during modeling

inferCNVAlleleBayesNet <- function(mcmc_allele,
                                   output_path,
                                   cores = 5){
  
  if(!dir.exists(file.path(output_path))){
    dir.create(file.path(output_path), recursive = T)
  }
  
  start_time <- Sys.time()
  if(class(mcmc_allele) == "MCMC_infercnv_allele_snp"){
    
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
  print(paste0("MCMC running time: ",
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
                                 thin = 1,
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
                                 thin = 1,
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
  
  if(ncol(theta_samples_df) == 3){
    colnames(theta_samples_df) <- c("del", "neutral","amp")
  } else if(ncol(theta_samples_df) == 7){
    colnames(theta_samples_df) <- c("del x2", "del x1",
                                    "neutral",
                                    "amp x1", "amp x2", "amp x3",
                                    "cnLOH")
  } else if(ncol(theta_samples_df) == 4){
    colnames(theta_samples_df) <- c("del", 
                                    "neutral",
                                    "amp", 
                                    "cnLOH")
    
  } else if(ncol(theta_samples_df) == 6){
    colnames(theta_samples_df) <- c("del x2", "del x1",
                                    "neutral",
                                    "amp x1", "amp x2", "amp x3")
  } else{colnames(theta_samples_df) <- c("del", "neutral")}
     
  theta_samples_df <- theta_samples_df %>% reshape::melt()
  colnames(theta_samples_df) <- c("index", "State", "Prob")
  if(length(unique(theta_samples_df$State)) == 3){
    theta_samples_df$State <- factor(theta_samples_df$State,
                                     levels = c("del", "neutral","amp"))
  } else if(length(unique(theta_samples_df$State)) == 7){
    theta_samples_df$State <- factor(theta_samples_df$State,
                                     levels = c("del x2", "del x1", 
                                                "neutral",
                                                "amp x1", "amp x2", "amp x3",
                                                "cnLOH"))
  } else if(length(unique(theta_samples_df$State)) == 4){
    theta_samples_df$State <- factor(theta_samples_df$State,
                                     levels = c("del", "neutral","amp","cnLOH"))
  } else if(length(unique(theta_samples_df$State)) == 6){
    theta_samples_df$State <- factor(theta_samples_df$State,
                                     levels = c("del x2", "del x1", 
                                                "neutral",
                                                "amp x1", "amp x2", "amp x3"))
  } else{
    theta_samples_df$State <- factor(theta_samples_df$State,
                                     levels = c("del", "neutral"))
  }
  
  stat_box_data <- function(y, upper_limit = max(theta_samples_df$Prob) * 1.15) {
    return( 
      data.frame(
        y = 0.95 * upper_limit,
        label = paste('count =', length(y), '\n',
                      'mean =', round(mean(y), 3), '\n')
      )
    )
  }
  
  pdf(file = file.path(output_path,
                       paste0(region,"_DiagnosticPlots.pdf")), onefile = TRUE)
  coda::autocorr.plot(coda::as.mcmc.list(theta_samples))
  coda::gelman.plot(coda::mcmc.list(theta_samples))
  dev.off()
  
  pdf(file = file.path(output_path,
                       paste0(region,"_DensityPlots.pdf")), onefile = TRUE)
  plot(samples)
  dev.off()
  
  prob_plot <- ggplot(theta_samples_df,
                      aes(y=Prob, x=State, color = State)) +
    geom_boxplot() +
    #ylim(c(0,1)) +
    stat_summary(fun.data = stat_box_data, 
                 geom = "text", 
                 hjust = 0.5,
                 vjust = 0.9) +
    theme_bw() + 
    ggtitle(region)
  # dev.off()
  ggsave(file.path(output_path,
                   paste0(region,"_CNVProbPlots.pdf")),
         prob_plot)
}

# inferCNVAlleleBayesNet <- function(infercnv_allele_obj,
#                                    infercnv_allele_hmm,
#                                    infercnv_allele_cellindex){
#   
#   #### testing codes for hmm inputs
#   genes.of.interest_list <- sapply(infercnv_allele_hmm, function(x)
#     unique(infercnv_allele_obj@SNP_info[x]$gene_name))
#   
#   ## associate each gene factor with a set of snps
#   genes2snps.dict_list <- lapply(seq_along(genes.of.interest_list), function(x) {
#     #browser()
#     genes2snps.dict <- lapply(seq_along(genes.of.interest_list[[x]]), function(y){
#       tmp_snp <- infercnv_allele_obj@SNP_info[infercnv_allele_hmm[[x]]]
#       tmp_name <- names(tmp_snp)
#       tmp_name <- tmp_name[which(tmp_snp$gene_name %in% genes.of.interest_list[[x]][y])]
#       })
#     names(genes2snps.dict) <- genes.of.interest_list[[x]]
#     return(genes2snps.dict)
#   })
#   ####
#   
#   ## initialize the arrays for processing
#   summary_array <- lapply(seq_along(infercnv_allele_hmm), function(x){
#     list("r.maf" = infercnv_allele_obj@count.data[infercnv_allele_hmm[[x]],
#                                                   infercnv_allele_cellindex[[x]]],
#          "n.sc" = infercnv_allele_obj@coverage.data[infercnv_allele_hmm[[x]],
#                                                     infercnv_allele_cellindex[[x]]],
#          "l.maf" = rowSums(infercnv_allele_obj@count.data[infercnv_allele_hmm[[x]],
#                                                           infercnv_allele_cellindex[[x]]] > 0),
#          "n.bulk" = rowSums(infercnv_allele_obj@coverage.data[infercnv_allele_hmm[[x]],
#                                                               infercnv_allele_cellindex[[x]]] > 0))
#   })
#   
#   ## input to MCMC
#   MCMC_inputs <- lapply(seq_along(summary_array), 
#                         function(x) populate_array(summary_array[[x]],
#                                                    genes2snps.dict_list[[x]]))
# 
#   #return(MCMC_inputs)
#   futile.logger::flog.info("Creating a new MCMC_allele_obj")
#   mcmc_allele_obj <- new(Class = "MCMC_infercnv_allele",
#                          bugs_model = system.file("BUGS_SNP_Model", package = "infercnv"),
#                          n.array = lapply(MCMC_inputs, function(x) x$input.array),
#                          n.meta = lapply(MCMC_inputs, function(x) x$input.metadata))
#   
#   validate_mcmc_infercnv_allele_obj(mcmc_allele_obj)
#   
#   ## mcmc run
#   futile.logger::flog.info("Start running MCMC")
#   start_time <- Sys.time()
#   mcmc_allele_obj <- runAlleleMCMC(mcmc_allele_obj)
#   end_time <- Sys.time()
#   futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
#   
#   return(mcmc_allele_obj)
# }

## internal function to get arrays necessary to the simulation                                    
# populate_array <- function(array_list, genes2snps.dict){
#   #browser()
#   I.j <- unlist(lapply(genes2snps.dict, length))
#   numGenes <- length(genes2snps.dict)
#   numSnpsPerGene <- max(I.j)
#   numCells <- ncol(array_list$r.maf)
#   
#   r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
#   for(i in seq_len(numGenes)) {
#     snpst <- genes2snps.dict[[i]]
#     for(s in seq_along(snpst)) {
#       r.array[i,s,] <- array_list$r.maf[snpst[s],]
#     }
#   }
#   n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
#   for(i in seq_len(numGenes)) {
#     snpst <- genes2snps.dict[[i]]
#     for(s in seq_along(snpst)) {
#       n.sc.array[i,s,] <- array_list$n.sc[snpst[s],]
#     }
#   }
#   l.array <- array(0, c(numGenes, numSnpsPerGene))
#   for(i in seq_len(numGenes)) {
#     snpst <- genes2snps.dict[[i]]
#     for(s in seq_along(snpst)) {
#       l.array[i,s] <- array_list$l.maf[snpst[s]]
#     }
#   }
#   n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
#   for(i in seq_len(numGenes)) {
#     snpst <- genes2snps.dict[[i]]
#     for(s in seq_along(snpst)) {
#       n.bulk.array[i,s] <- array_list$n.bulk[snpst[s]]
#     }
#   }
#   
#   return(list(input.array = list("r.array" = r.array,
#                                  "n.sc.array" = n.sc.array,
#                                  "l.array" = l.array,
#                                  "n.bulk.array" = n.bulk.array),
#               input.metadata = list("nGenes" = numGenes,
#                                     "nCells" = numCells,
#                                     "nsnps2genes" = I.j)))
# }
# 
# ## validation function for creating mcmc_infercnv_allele obj
# validate_mcmc_infercnv_allele_obj <- function(mcmc_allele_obj){
#   flog.info("validating mcmc_infercnv_allele_obj")
#   
#   if(isTRUE(length(mcmc_allele_obj@n.array) == length(mcmc_allele_obj@n.meta))){
#     return()
#   } else{
#     flog.error("hmm.... length(mcmc_allele_obj@n.array != 
#                rownames(mcmc_allele_obj@n.meta))")
#   }
#   
#   broken.mcmc_infercnv_allele_obj = mcmc_allele_obj
#   save(broken.mcmc_infercnv_allele_obj, file="broken.mcmc_infercnv_allele_obj")
#   stop("Problem detected w/ mcmc_infercnv_allele_obj")
# }
# 
# 
# ## internal function to run allele-based only mcmc
# runAlleleMCMC <- function(mcmc_allele_obj, pe = 0.1, mono = 0.7, n.iter = 1e3){
#   
#   mcmc <- mclapply(seq_along(mcmc_allele_obj@n.array), function(x) {
#     
#     data <- list(
#       'l' = mcmc_allele_obj@n.array[[x]]$l.array,
#       'r' = mcmc_allele_obj@n.array[[x]]$r.array,
#       'n.bulk' = mcmc_allele_obj@n.array[[x]]$n.bulk.array,
#       'n.sc' = mcmc_allele_obj@n.array[[x]]$n.sc.array,
#       'J' = mcmc_allele_obj@n.meta[[x]]$nGenes,  # how many genes
#       'K' = mcmc_allele_obj@n.meta[[x]]$nCells,  # how many cells
#       'I.j' = mcmc_allele_obj@n.meta[[x]]$nsnps2genes,
#       'pseudo' = pe,
#       'mono' = mono)
#     
#     model <- rjags::jags.model(mcmc_allele_obj@bugs_model, 
#                                data=data, n.chains=4, n.adapt=300)
#     update(model, 300)
#     cat('Done modeling!')
#     
#     parameters <- c('alpha', 'S')
#     samples <- coda.samples(model, parameters, n.iter=n.iter)
#   }, mc.cores = 6)
#   
#   mcmc_allele_obj@posterior_prob <- mcmc
#   
#   return(mcmc_allele_obj)
# }
# 
# ## internal function to run combined-based mcmc
# ## gr disjon-based grange obj based on allele/gene boundaries
# ## gene_annot gene annotation from infercnv gene-based obj
# runCombinedMCMC <- function(gr, gene_annot,
#                             infercnv_allele_obj, mcmc_gene_obj,
#                             mono = 0.7, pe = 0.1, n.iter = 1e3,
#                             mode = mode){
#   
#   lapply(sapply(gr$region, list), function(x) {
#     #browser()
#     snp_map_index <- infercnv_allele_obj@SNP_info %over% gr[gr$region == x]
#     gene_map_index <- gene_annot %over% gr[gr$region == x]
#     
#     if(sum(snp_map_index) <= 4 | sum(gene_map_index) <= 2){
#       return()
#     }
#     
#     r.maf <- infercnv_allele_obj@count.data[snp_map_index, 
#                                             unlist(infercnv_allele_obj@observation_grouped_cell_indices),
#                                             drop = F]
#     n.sc = infercnv_allele_obj@coverage.data[snp_map_index,
#                                              unlist(infercnv_allele_obj@observation_grouped_cell_indices),
#                                              drop = F]
#     
#     if(mode %in% c("combined", "allele", "HB_allele", "gene")){
#       
#       l.maf = rowSums(r.maf > 0)
#       n.bulk = rowSums(n.sc > 0)
#       
#       geneFactor <- infercnv_allele_obj@SNP_info$gene_name[snp_map_index]
#       genes.of.interest <- unique(geneFactor)
#       genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
#         names(infercnv_allele_obj@SNP_info[snp_map_index])[which(geneFactor %in% genes.of.interest[i])]
#       })
#       names(genes2snps.dict) <- genes.of.interest
#       
#       ## Convert to multi-dimensions based on j
#       I.j <- unlist(lapply(genes2snps.dict, length))
#       numGenes <- length(genes2snps.dict)
#       numSnpsPerGene <- max(I.j)
#       numCells <- ncol(r.maf)## original name error --Rongting
#       ## j, i, k
#       r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
#       for(i in seq_len(numGenes)) {
#         snpst <- genes2snps.dict[[i]]
#         for(s in seq_along(snpst)) {
#           r.array[i,s,] <- r.maf[snpst[s],]
#         }
#       }
#       n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
#       for(i in seq_len(numGenes)) {
#         snpst <- genes2snps.dict[[i]]
#         for(s in seq_along(snpst)) {
#           n.sc.array[i,s,] <- n.sc[snpst[s],]
#         }
#       }
#       l.array <- array(0, c(numGenes, numSnpsPerGene))
#       for(i in seq_len(numGenes)) {
#         snpst <- genes2snps.dict[[i]]
#         for(s in seq_along(snpst)) {
#           l.array[i,s] <- l.maf[snpst[s]]
#         }
#       }
#       n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
#       for(i in seq_len(numGenes)) {
#         snpst <- genes2snps.dict[[i]]
#         for(s in seq_along(snpst)) {
#           n.bulk.array[i,s] <- n.bulk[snpst[s]]
#         }
#       }
#       
#       gexp <- mcmc_gene_obj@expr.data[gene_map_index,
#                                       unlist(infercnv_allele_obj@observation_grouped_cell_indices),
#                                       drop = F]
#       
#       data <- list(
#         'l' = l.array,
#         'r' = r.array,
#         'n.bulk' = n.bulk.array,
#         'n.sc' = n.sc.array,
#         'J' = numGenes,  # how many genes
#         'K' = numCells,  # how many cells
#         'I.j' = I.j,
#         'pseudo' = pe,
#         'mono' = mono,
#         'gexp' = gexp,
#         'JJ' = nrow(gexp),
#         "mu" = mcmc_gene_obj@mu,
#         "sig" = mcmc_gene_obj@sig)
#       
#     } else if(mode == "snp2gene"){
#       
#       test_gene.rmaf <- lapply(seq_len(ncol(r.maf)), 
#                                function(x){
#                                  tapply(r.maf[,x], 
#                                         infercnv_allele_obj@SNP_info$gene_name[snp_map_index], 
#                                         mean)
#                                }) 
#       test_gene.rmaf <- do.call(cbind, test_gene.rmaf)
#       colnames(test_gene.rmaf) <- colnames(r.maf)
#       r.array <- round(test_gene.rmaf)
#       
#       test_gene.n.sc <- lapply(seq_len(ncol(n.sc)), 
#                                function(x){
#                                  tapply(n.sc[,x], 
#                                         infercnv_allele_obj@SNP_info$gene_name[snp_map_index], 
#                                         mean)
#                                }) 
#       test_gene.n.sc <- do.call(cbind, test_gene.n.sc)
#       colnames(test_gene.n.sc) <- colnames(n.sc)
#       n.sc.array <- round(test_gene.n.sc)
#       
#       l.array <- rowSums(r.array > 0)
#       n.bulk.array <- rowSums(n.sc.array > 0)
#       
#       data <- list(
#         'l' = l.array,
#         'r' = r.array,
#         'n.bulk' = n.bulk.array,
#         'n.sc' = n.sc.array,
#         'J' = nrow(r.array),  # how many genes
#         'K' = ncol(r.array),  # how many cells
#         'pseudo' = pe,
#         'mono' = mono)
#       
#     }
#     
#     if(mode %in% c("allele","snp2gene")){
#       inits <- list(
#         list(epsilon = rep(1, data$K)),
#         list(epsilon = rep(1, data$K)),
#         list(epsilon = rep(2, data$K)),
#         list(epsilon = rep(2, data$K))
#       )
#       n.chain <- 4
#     } else{
#       inits <- list(
#         list(epsilon = rep(1, data$K)),
#         list(epsilon = rep(2, data$K)),
#         list(epsilon = rep(3, data$K))
#       )
#       n.chain <- 3
#     }
#     
#     if (mode == "HB_allele"){
#       
#       model <- rjags::jags.model(system.file("BUGS_SNP_Model", package = "infercnv"), 
#                                  data=data, n.chains=3, n.adapt=500)
#       update(model, 200)
#       parameters <- c('alpha', 'S')
#       samples <- coda.samples(model, parameters, n.iter=n.iter)
#       
#     } else{
#       
#       model <- rjags::jags.model(ifelse(mode == "combined",
#                                         system.file("BUGS_Combined_Model_i3", package = "infercnv"),
#                                         ifelse(mode == "gene", 
#                                                system.file("BUGS_Mixture_Model_i3_gene_test", package = "infercnv"),
#                                                ifelse(mode == "allele", 
#                                                       system.file("BUGS_SNP_Model_i2_test", package = "infercnv"),
#                                                       system.file("BUGS_SNP2Gene_Model", package = "infercnv")))),
#                                  data = data,
#                                  inits = inits,
#                                  n.chains = n.chain, 
#                                  n.adapt = 500)
#       
#       update(model, 200)
#       parameters <- c('theta', 'epsilon')
#       samples <- coda.samples(model, parameters, n.iter=n.iter)
#       
#     }
#     return(samples)
#   })
# }
# 
# ## the main function to run combined method -- under development
# ## infercnv_allele_obj
# ## infercnv_allele_hmm a list of index of HMM boundaries based on the allele info
# ## mcmc_gene_obj gene based mcmc obj (from maxwell's work)
# 
# inferCNVCombinedBayesNet <- function(infercnv_allele_obj,
#                                      infercnv_allele_hmm, 
#                                      mcmc_gene_obj,
#                                      #combine_method = c("disjoin", "intersect"),
#                                      mono = 0.7, pe = 0.1, n.iter = 1e3,
#                                      cores = 6, mode = c("combined","gene", 
#                                                          "allele","HB_allele",
#                                                          "snp2gene")){
#   mode = match.arg(mode)
#   
#   snp_gr <- do.call("c", lapply(infercnv_allele_hmm, function(x) 
#     range(infercnv_allele_obj@SNP_info[x])))
#   
#   gene_annot <- GRanges(seqnames = mcmc_gene_obj@gene_order[[C_CHR]],
#                         IRanges(mcmc_gene_obj@gene_order[[C_START]], 
#                                 mcmc_gene_obj@gene_order[[C_STOP]]))
#   gene_index <- c()
#   gene_index <- sapply(mcmc_gene_obj@cell_gene, 
#                        function(x) gene_index <- c(gene_index, x$Genes))
#   gene_gr <- do.call("c", lapply(gene_index, function(x) 
#     range(gene_annot[x])))
#   
#   union_gr <- c(snp_gr, gene_gr) %>% GenomicRanges::reduce()
#   # union_gr <- c(snp_gr, gene_gr) %>% disjoin()
#   # union_gr$region <- ifelse(union_gr %over% snp_gr, 
#   #                           ifelse(union_gr %over% gene_gr,
#   #                                  "common_region",
#   #                                  ifelse(overlapsAny(union_gr, snp_gr, type = "equal"),
#   #                                         "snp_only","disjoin_snp")),
#   #                           ifelse(overlapsAny(union_gr, gene_gr, type = "equal"),
#   #                                  "gene_only", "disjoin_gene"))
#   # union_gr$region <- paste0(union_gr$region, "_", seq_along(union_gr$region))
#   union_gr$region <- paste0("region_", seq_along(union_gr))
#   union_gr <- union_gr %>% split(seqnames(union_gr), 
#                                  drop = T) %>% GRangesList()
#   
#   futile.logger::flog.info("Start running MCMC")
#   futile.logger::flog.info(paste0("Running mode: ", mode))
#   start_time <- Sys.time()
#   
#   sim_res <- mclapply(union_gr, function(x)
#     #browser()
#     runCombinedMCMC(x, gene_annot,
#                     infercnv_allele_obj, mcmc_gene_obj,
#                     mono = 0.7, pe = 0.1, n.iter = 1e3,
#                     mode = mode),
#     mc.cores = cores)
#   
#   end_time <- Sys.time()
#   futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
#   
#   return(sim_res)
# }
