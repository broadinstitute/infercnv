#' MCMC infercnv_allele class -- under development
#' 
#' @description This class extends the functionality of infercnv_allele class
#' aiming to coorporate MCMC object leveraging Bayesian model
#' 
#' Slots in the MCMC_infercnv_allele object (besides infercnv) include:
#' 
#' @slot bugs_model a string path to where the bug file locates
#' 
#' @slot n.array the list of arrays containing allele count/ ref count etc.
#' 
#' @slot n.meta metadata info
#' 
#' @slot posterior_prob the list obj containing posterior prob from bayesian model
#' 
#' @export
#' 
MCMC_infercnv_allele <- methods::setClass(
  "MCMC_infercnv_allele", 
  slots = c(bugs_model = "character",
            n.array = "list",
            n.meta = "list",
            posterior_prob = "list"),
  contains = "infercnv_allele")

#' @title inferCNVAlleleBayesNet: Run Bayesian NetworkModel To Obtain Posterior Probabilities For HMM Predicted States
#' @description Run MCMC method to get posterior probability of cells having Deletions/LOHs
#' @param infercnv_allele_obj an infercnv_allele obj
#' @param infercnv_allele_hmm a list of index of HMM boundaries
#' @param infercnv_allele_cellindex an index of cell annotation

inferCNVAlleleBayesNet <- function(infercnv_allele_obj,
                                   infercnv_allele_hmm,
                                   infercnv_allele_cellindex){
  
  #### testing codes for hmm inputs
  genes.of.interest_list <- sapply(infercnv_allele_hmm, function(x)
    unique(infercnv_allele_obj@SNP_info[x]$gene_name))
  
  ## associate each gene factor with a set of snps
  genes2snps.dict_list <- lapply(seq_along(genes.of.interest_list), function(x) {
    #browser()
    genes2snps.dict <- lapply(seq_along(genes.of.interest_list[[x]]), function(y){
      tmp_snp <- infercnv_allele_obj@SNP_info[infercnv_allele_hmm[[x]]]
      tmp_name <- names(tmp_snp)
      tmp_name <- tmp_name[which(tmp_snp$gene_name %in% genes.of.interest_list[[x]][y])]
      })
    names(genes2snps.dict) <- genes.of.interest_list[[x]]
    return(genes2snps.dict)
  })
  ####
  
  ## initialize the arrays for processing
  summary_array <- lapply(seq_along(infercnv_allele_hmm), function(x){
    list("r.maf" = infercnv_allele_obj@count.data[infercnv_allele_hmm[[x]],
                                                  infercnv_allele_cellindex[[x]]],
         "n.sc" = infercnv_allele_obj@coverage.data[infercnv_allele_hmm[[x]],
                                                    infercnv_allele_cellindex[[x]]],
         "l.maf" = rowSums(infercnv_allele_obj@count.data[infercnv_allele_hmm[[x]],
                                                          infercnv_allele_cellindex[[x]]] > 0),
         "n.bulk" = rowSums(infercnv_allele_obj@coverage.data[infercnv_allele_hmm[[x]],
                                                              infercnv_allele_cellindex[[x]]] > 0))
  })
  
  ## input to MCMC
  MCMC_inputs <- lapply(seq_along(summary_array), 
                        function(x) populate_array(summary_array[[x]],
                                                   genes2snps.dict_list[[x]]))

  #return(MCMC_inputs)
  futile.logger::flog.info("Creating a new MCMC_allele_obj")
  mcmc_allele_obj <- new(Class = "MCMC_infercnv_allele",
                         bugs_model = system.file("BUGS_SNP_Model", package = "infercnv"),
                         n.array = lapply(MCMC_inputs, function(x) x$input.array),
                         n.meta = lapply(MCMC_inputs, function(x) x$input.metadata))
  
  validate_mcmc_infercnv_allele_obj(mcmc_allele_obj)
  
  ## mcmc run
  futile.logger::flog.info("Start running MCMC")
  start_time <- Sys.time()
  mcmc_allele_obj <- runAlleleMCMC(mcmc_allele_obj)
  end_time <- Sys.time()
  futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
  
  return(mcmc_allele_obj)
}

## internal function to get arrays necessary to the simulation                                    
populate_array <- function(array_list, genes2snps.dict){
  #browser()
  I.j <- unlist(lapply(genes2snps.dict, length))
  numGenes <- length(genes2snps.dict)
  numSnpsPerGene <- max(I.j)
  numCells <- ncol(array_list$r.maf)
  
  r.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      r.array[i,s,] <- array_list$r.maf[snpst[s],]
    }
  }
  n.sc.array <- array(0, c(numGenes, numSnpsPerGene, numCells))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      n.sc.array[i,s,] <- array_list$n.sc[snpst[s],]
    }
  }
  l.array <- array(0, c(numGenes, numSnpsPerGene))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      l.array[i,s] <- array_list$l.maf[snpst[s]]
    }
  }
  n.bulk.array <- array(0, c(numGenes, numSnpsPerGene))
  for(i in seq_len(numGenes)) {
    snpst <- genes2snps.dict[[i]]
    for(s in seq_along(snpst)) {
      n.bulk.array[i,s] <- array_list$n.bulk[snpst[s]]
    }
  }
  
  return(list(input.array = list("r.array" = r.array,
                                 "n.sc.array" = n.sc.array,
                                 "l.array" = l.array,
                                 "n.bulk.array" = n.bulk.array),
              input.metadata = list("nGenes" = numGenes,
                                    "nCells" = numCells,
                                    "nsnps2genes" = I.j)))
}

## validation function for creating mcmc_infercnv_allele obj
validate_mcmc_infercnv_allele_obj <- function(mcmc_allele_obj){
  flog.info("validating mcmc_infercnv_allele_obj")
  
  if(isTRUE(length(mcmc_allele_obj@n.array) == length(mcmc_allele_obj@n.meta))){
    return()
  } else{
    flog.error("hmm.... length(mcmc_allele_obj@n.array != 
               rownames(mcmc_allele_obj@n.meta))")
  }
  
  broken.mcmc_infercnv_allele_obj = mcmc_allele_obj
  save(broken.mcmc_infercnv_allele_obj, file="broken.mcmc_infercnv_allele_obj")
  stop("Problem detected w/ mcmc_infercnv_allele_obj")
}


## internal function to run allele-based only mcmc
runAlleleMCMC <- function(mcmc_allele_obj, pe = 0.1, mono = 0.7, n.iter = 1e3){
  
  mcmc <- mclapply(seq_along(mcmc_allele_obj@n.array), function(x) {
    
    data <- list(
      'l' = mcmc_allele_obj@n.array[[x]]$l.array,
      'r' = mcmc_allele_obj@n.array[[x]]$r.array,
      'n.bulk' = mcmc_allele_obj@n.array[[x]]$n.bulk.array,
      'n.sc' = mcmc_allele_obj@n.array[[x]]$n.sc.array,
      'J' = mcmc_allele_obj@n.meta[[x]]$nGenes,  # how many genes
      'K' = mcmc_allele_obj@n.meta[[x]]$nCells,  # how many cells
      'I.j' = mcmc_allele_obj@n.meta[[x]]$nsnps2genes,
      'pseudo' = pe,
      'mono' = mono)
    
    model <- rjags::jags.model(mcmc_allele_obj@bugs_model, 
                               data=data, n.chains=4, n.adapt=300)
    update(model, 300)
    cat('Done modeling!')
    
    parameters <- c('alpha', 'S')
    samples <- coda.samples(model, parameters, n.iter=n.iter)
  }, mc.cores = 6)
  
  mcmc_allele_obj@posterior_prob <- mcmc
  
  return(mcmc_allele_obj)
}

## internal function to run combined-based mcmc
## gr disjon-based grange obj based on allele/gene boundaries
## gene_annot gene annotation from infercnv gene-based obj
runCombinedMCMC <- function(gr, gene_annot,
                            infercnv_allele_obj, mcmc_gene_obj,
                            mono = 0.7, pe = 0.1, n.iter = 1e3){
  
  lapply(seq_along(gr), function(x) {
    #browser()
    snp_map_index <- infercnv_allele_obj@SNP_info %over% gr[x]
    gene_map_index <- gene_annot %over% gr[x]
    
    if(mean(snp_map_index) == 0 | mean(gene_map_index) == 0){
      return()
    }
    
    r.maf <- infercnv_allele_obj@count.data[snp_map_index, 
                                            unlist(infercnv_allele_obj@observation_grouped_cell_indices),
                                            drop = F]
    n.sc = infercnv_allele_obj@coverage.data[snp_map_index,
                                             unlist(infercnv_allele_obj@observation_grouped_cell_indices),
                                             drop = F]
    l.maf = rowSums(r.maf > 0)
    n.bulk = rowSums(n.sc > 0)
    
    geneFactor <- infercnv_allele_obj@SNP_info$gene_name[snp_map_index]
    genes.of.interest <- unique(geneFactor)
    genes2snps.dict <- lapply(seq_along(genes.of.interest), function(i) {
      names(infercnv_allele_obj@SNP_info[snp_map_index])[which(geneFactor %in% genes.of.interest[i])]
    })
    names(genes2snps.dict) <- genes.of.interest
    
    ## Convert to multi-dimensions based on j
    I.j <- unlist(lapply(genes2snps.dict, length))
    numGenes <- length(genes2snps.dict)
    numSnpsPerGene <- max(I.j)
    numCells <- ncol(r.maf)## original name error --Rongting
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
    
    gexp <- mcmc_gene_obj@expr.data[gene_map_index,
                                      unlist(infercnv_allele_obj@observation_grouped_cell_indices),
                                      drop = F]
    
    data <- list(
      'l' = l.array,
      'r' = r.array,
      'n.bulk' = n.bulk.array,
      'n.sc' = n.sc.array,
      'J' = numGenes,  # how many genes
      'K' = numCells,  # how many cells
      'I.j' = I.j,
      'pseudo' = pe,
      'mono' = mono,
      'gexp' = gexp,
      'JJ' = nrow(gexp),
      "mu" = mcmc_gene_obj@mu,
      "sig" = mcmc_gene_obj@sig)
    
    model <- rjags::jags.model(system.file("BUGS_Combined_Model_i3", package = "infercnv"), 
                               data=data, n.chains=4, n.adapt=300)
    update(model, 300)
    parameters <- c('theta', 'epsilon')
    samples <- coda.samples(model, parameters, n.iter=n.iter)
    
    return(samples)
  })
}

## the main function to run combined method -- under development
## infercnv_allele_obj
## infercnv_allele_hmm a list of index of HMM boundaries based on the allele info
## mcmc_gene_obj gene based mcmc obj (from maxwell's work)

inferCNVCombinedBayesNet <- function(infercnv_allele_obj,
                                     infercnv_allele_hmm, 
                                     mcmc_gene_obj,
                                     #combine_method = c("disjoin", "intersect"),
                                     mono = 0.7, pe = 0.1, n.iter = 1e3,
                                     cores = 6){
  
  snp_gr <- do.call("c", lapply(infercnv_allele_hmm, function(x) 
    range(infercnv_allele_obj@SNP_info[x])))
  
  gene_annot <- GRanges(seqnames = mcmc_gene_obj@gene_order[[C_CHR]],
                        IRanges(mcmc_gene_obj@gene_order[[C_START]], 
                                mcmc_gene_obj@gene_order[[C_STOP]]))
  gene_index <- c()
  gene_index <- sapply(mcmc_gene_obj@cell_gene, 
                       function(x) gene_index <- c(gene_index, x$Genes))
  gene_gr <- do.call("c", lapply(gene_index, function(x) 
    range(gene_annot[x])))
  
  union_gr <- c(snp_gr, gene_gr)
  union_gr <- union_gr %>% split(seqnames(union_gr), 
                                 drop = T) %>% GRangesList() %>% disjoin()
  
  futile.logger::flog.info("Start running MCMC")
  start_time <- Sys.time()
  
  sim_res <- mclapply(union_gr, function(x)
    #browser()
    runCombinedMCMC(x, gene_annot,
                    infercnv_allele_obj, mcmc_gene_obj,
                    mono = 0.7, pe = 0.1, n.iter = 1e3),
    mc.cores = cores)
  
  end_time <- Sys.time()
  futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
  
  return(sim_res)
}
