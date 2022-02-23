#' MCMC infercnv_allele class
#' 
#' @description This class extends the functionality of infercnv_allele class
#' aiming to coorporate MCMC object leveraging Bayesian model
#' 
#' Slots in the MCMC_infercnv_allele object (besides infercnv) include:
#' 
#' @slot bugs_model
#' 
#' @slot n.array
#' 
#' @slot n.meta
#' 
#' @slot alpha_prob
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

inferCNVAlleleBayesNet <- function(infercnv_allele_obj,
                                   infercnv_allele_hmm){
  
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
  summary_array <- lapply(infercnv_allele_hmm, function(x){
    list("r.maf" = infercnv_allele_obj@count.data[x,],
         "n.sc" = infercnv_allele_obj@coverage.data[x,],
         "l.maf" = rowSums(infercnv_allele_obj@count.data[x,] > 0),
         "n.bulk" = rowSums(infercnv_allele_obj@coverage.data[x,] > 0))
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