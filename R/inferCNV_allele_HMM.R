#' @title allele_HMM_predict_CNV_via_HMM_on_tumor_subclusters
#'
#' @description use the allele HMM for predicting CNV at the level of tumor subclusters
#'
#' @return infercnv_allele_obj where infercnv_allele_obj@expr.data contains state assignments.
#'
#' @keywords internal
#' 
#' @noRd
allele_HMM_predict_CNV_via_HMM_on_tumor_subclusters <- function(infercnv_allele_obj,
                                                                t = 1e-6, pd = 0.1, pn = 0.45,
                                                                min.num.snps = 5, trim = 0.1){
  ## pre-check for allele data
  if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
    flog.info("Initializing the lesser allele fraction ...")
    infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
  }
  
  flog.info("predict_allele_CNV_via_HMM_on_tumor_subclusters")
  
  chrs = unique(infercnv_allele_obj@gene_order$chr)
  gene_order = infercnv_allele_obj@gene_order
  lesser.data <- infercnv_allele_obj@count.data
  coverage.data <- infercnv_allele_obj@coverage.data
  
  ## initialize hmm states for allele data
  hmm.allele.data <- matrix(0,
                            nrow = nrow(lesser.data),
                            ncol = ncol(lesser.data))
  rownames(hmm.allele.data) <- rownames(lesser.data)
  colnames(hmm.allele.data) <- colnames(lesser.data)
  
  tumor_subclusters <- unlist(infercnv_allele_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
  
  HMM_output <- c()
  cell_index <- c()
  ## add the normals, so they get predictions too:
  #tumor_subclusters <- c(tumor_subclusters, infercnv_allele_obj@reference_grouped_cell_indices)
  
  ## run HMM across chromosomes
  lapply(chrs, function(chr){
    
    chr_snp_idx = which(gene_order$chr == chr)
    
    lapply(tumor_subclusters, function(tumor_subcluster_cells_idx){
      
      if (length(tumor_subcluster_cells_idx) > 1){
        
        if(mean(tumor_subcluster_cells_idx %in% 
                unlist(infercnv_allele_obj@observation_grouped_cell_indices)) == 1){
          status <- "tumor"
        }
        if(mean(tumor_subcluster_cells_idx %in% 
                unlist(infercnv_allele_obj@reference_grouped_cell_indices)) == 1){
          status <- "normal"
        }
        
        lesser.data <- lesser.data[chr_snp_idx, tumor_subcluster_cells_idx,drop=FALSE]
        coverage.data <- coverage.data[chr_snp_idx, tumor_subcluster_cells_idx,drop=FALSE]
        
        mafl <- rowSums(lesser.data > 0)
        sizel <- rowSums(coverage.data > 0)
        
        ## change point
        delta <- c(0.5, 0.5)
        z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                byrow=TRUE, nrow=2), 
                   delta, "binom", list(prob=c(pd, pn)), 
                   list(size=sizel), discrete=TRUE)
        #results <- Viterbi(z)
        results <- Viterbi.dthmm.allele.adj(z)
        boundsnps <- rownames(lesser.data)[results == 1]
        
        ## vote
        vote <- rep(0, nrow(lesser.data))
        names(vote) <- rownames(lesser.data)
        vote[boundsnps] <- 1
        
        if(max(vote) == 0) {
          flog.info(sprintf('Exiting; no new bound SNPs found at %s in %s cells ...', 
                            chr, status))
          return() ## exit iteration, no more bound SNPs found
        }
        
        vote[vote > 0] <- 1
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.snps.cont <- rep(0, length(vote))
        names(bound.snps.cont) <- names(vote)
        
        for(i in 2:length(vote)) {
          if(vote[i] >= mv & vote[i] == vote[i-1]) {
            bound.snps.cont[i] <- cs
          } else {
            cs <- cs + 1
          }
        }
        
        tb <- table(bound.snps.cont)
        tbv <- as.vector(tb)
        names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0
        
        ## all detected deletions have fewer than 5 SNPs...reached the end
        tbv <- tbv[tbv >= min.num.snps]
        if(length(tbv)==0) {
          flog.info(sprintf('Exiting; less than %s new bound SNPs found at %s in %s cells ...', 
                            min.num.snps, chr, status))
          return()
        }
        
        HMM_info <- lapply(names(tbv), function(ti) {
          
          bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
          
          ## trim
          bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
          cell_index <<- c(cell_index, list(tumor_subcluster_cells_idx))
          
          return(bound.snps.new)
          
        })
        HMM_output <<- c(HMM_output, HMM_info)
        HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
        
        flog.info(sprintf("Done extracting HMM regions at %s in %s cells ...", 
                          chr, status))
        print(HMM_region)
        
        if(!is.null(HMM_region)){
          hmm.allele.data[infercnv_allele_obj@SNP_info %over% HMM_region,
                          tumor_subcluster_cells_idx] <<- -1
        }
      }
    })
  })
  infercnv_allele_obj@expr.data <- hmm.allele.data
  # return(list(infercnv_allele_obj = infercnv_allele_obj, 
  #             HMM_output = HMM_output,
  #             cell_index = cell_index))
  return(infercnv_allele_obj)
}

#' @title allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples
#'
#' @description use the allele HMM for predicting CNV at the level of whole tumor samples
#'
#' @return infercnv_allele_obj where infercnv_allele_obj@expr.data contains state assignments.
#'
#' @keywords internal
#' 
#' @noRd
allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples <- function(infercnv_allele_obj,
                                                                  t = 1e-6, pd = 0.1, pn = 0.45,
                                                                  min.num.snps = 5, trim = 0.1){
  ## pre-check for allele data
  if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
    flog.info("Initializing the lesser allele fraction ...")
    infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
  }
  
  flog.info("predict_allele_CNV_via_HMM_on_whole_tumor_samples")
  
  chrs = unique(infercnv_allele_obj@gene_order$chr)
  gene_order = infercnv_allele_obj@gene_order
  lesser.data <- infercnv_allele_obj@count.data
  coverage.data <- infercnv_allele_obj@coverage.data
  
  ## initialize hmm allele data
  hmm.allele.data <- matrix(0,
                            nrow = nrow(lesser.data),
                            ncol = ncol(lesser.data))
  rownames(hmm.allele.data) <- rownames(lesser.data)
  colnames(hmm.allele.data) <- colnames(lesser.data)
  
  ## add the normals, so they get predictions too:
  tumor_samples <- c(infercnv_allele_obj@observation_grouped_cell_indices, 
                     infercnv_allele_obj@reference_grouped_cell_indices)
  HMM_output <- c()
  cell_index <- c()
  
  ## run HMM across chromosomes
  lapply(chrs, function(chr){
    
    chr_snp_idx = which(gene_order$chr == chr)
    
    lapply(tumor_samples, function(tumor_sample_cells_idx){

      if(length(tumor_sample_cells_idx) > 1){
        
        if(mean(tumor_sample_cells_idx %in% 
                unlist(infercnv_allele_obj@observation_grouped_cell_indices)) == 1){
          status <- "tumor"
        }
        if(mean(tumor_sample_cells_idx %in% 
                unlist(infercnv_allele_obj@reference_grouped_cell_indices)) == 1){
          status <- "normal"
        }
        
        lesser.data <- lesser.data[chr_snp_idx, tumor_sample_cells_idx,drop=FALSE]
        coverage.data <- coverage.data[chr_snp_idx, tumor_sample_cells_idx,drop=FALSE]
        
        mafl <- rowSums(lesser.data > 0)
        sizel <- rowSums(coverage.data > 0)
        
        ## change point
        delta <- c(0.5, 0.5)
        z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                byrow=TRUE, nrow=2), 
                   delta, "binom", list(prob=c(pd, pn)), 
                   list(size=sizel), discrete=TRUE)
        #results <- Viterbi(z)
        results <- Viterbi.dthmm.allele.adj(z)
        boundsnps <- rownames(lesser.data)[results == 1]
        
        ## vote
        vote <- rep(0, nrow(lesser.data))
        names(vote) <- rownames(lesser.data)
        
        vote[boundsnps] <- 1
        
        if(max(vote) == 0) {
          flog.info(sprintf('Exiting; no new bound SNPs found at %s in %s cells ...', 
                            chr, status))
          return() ## exit iteration, no more bound SNPs found
        }
        
        vote[vote > 0] <- 1
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.snps.cont <- rep(0, length(vote))
        names(bound.snps.cont) <- names(vote)
        
        for(i in 2:length(vote)) {
          if(vote[i] >= mv & vote[i] == vote[i-1]) {
            bound.snps.cont[i] <- cs
          } else {
            cs <- cs + 1
          }
        }
        
        tb <- table(bound.snps.cont)
        tbv <- as.vector(tb)
        names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0
        
        ## all detected deletions have fewer than 5 SNPs...reached the end
        tbv <- tbv[tbv >= min.num.snps]
        if(length(tbv)==0) {
          flog.info(sprintf('Exiting; less than %s new bound SNPs found at %s in %s cells ...', 
                            min.num.snps, chr, status))
          return()
        }
        
        HMM_info <- lapply(names(tbv), function(ti) {
          
          bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
          
          ## trim
          bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
          
          cell_index <<- c(cell_index, list(tumor_sample_cells_idx))
          return(bound.snps.new)
        })
        HMM_output <<- c(HMM_output, HMM_info)
        HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
        
        flog.info(sprintf("Done extracting HMM regions at %s in %s cells ...", 
                          chr, status))
        print(HMM_region)
        
        if(!is.null(HMM_region)){
          hmm.allele.data[infercnv_allele_obj@SNP_info %over% HMM_region,
                          tumor_sample_cells_idx] <<- -1
        }
      }
    })
  })
  infercnv_allele_obj@expr.data <- hmm.allele.data
  # return(list(infercnv_allele_obj = infercnv_allele_obj, 
  #             HMM_output = HMM_output,
  #             cell_index = cell_index))
  return(infercnv_allele_obj)
}


################## deprecated 
allele_HMM_predict_CNV_via_HMM_on_tumor_subclusters_mod <- function(infercnv_allele_obj,
                                                                    t = 1e-6, pd = 0.1, pn = 0.45,
                                                                    min.num.snps = 5, trim = 0.1,
                                                                    min.traverse = 3){
  
  ## pre-check for allele data
  if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
    flog.info("Initializing the lesser allele fraction ...")
    infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
  }
  
  flog.info("predict_allele_CNV_via_HMM_on_tumor_subclusters")
  
  chrs = unique(infercnv_allele_obj@gene_order$chr)
  gene_order = infercnv_allele_obj@gene_order
  lesser.frac.data <- infercnv_allele_obj@expr.data
  lesser.data <- infercnv_allele_obj@count.data
  coverage.data <- infercnv_allele_obj@coverage.data
  
  ## initialize hmm states for allele data
  hmm.allele.data <- matrix(0,
                            nrow = nrow(lesser.data),
                            ncol = ncol(lesser.data))
  rownames(hmm.allele.data) <- rownames(lesser.data)
  colnames(hmm.allele.data) <- colnames(lesser.data)
  
  tumor_subclusters <- unlist(infercnv_allele_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
  
  HMM_output <- c()
  cell_index <- c()
  ## add the normals, so they get predictions too:
  #tumor_subclusters <- c(tumor_subclusters, infercnv_allele_obj@reference_grouped_cell_indices)
  
  ## run HMM across chromosomes
  lapply(chrs, function(chr){
    
    chr_snp_idx = which(gene_order$chr == chr)
    
    lapply(tumor_subclusters, function(tumor_subcluster_cells_idx){
      
      if (length(tumor_subcluster_cells_idx) > 1){
        
        if(mean(tumor_subcluster_cells_idx %in% 
                unlist(infercnv_allele_obj@observation_grouped_cell_indices)) == 1){
          status <- "tumor"
        }
        if(mean(tumor_subcluster_cells_idx %in% 
                unlist(infercnv_allele_obj@reference_grouped_cell_indices)) == 1){
          status <- "normal"
        }
        
        lesser.frac.data <- lesser.frac.data[chr_snp_idx,
                                             tumor_subcluster_cells_idx,drop=FALSE]
        lesser.data <- lesser.data[chr_snp_idx,
                                   tumor_subcluster_cells_idx,drop=FALSE]
        coverage.data <- coverage.data[chr_snp_idx,
                                       tumor_subcluster_cells_idx,drop=FALSE]
        
        d <- parallelDist(t(lesser.frac.data), 
                          method = "euclidean", threads = 5) # parallel dist
        hc <- hclust(d, method="ward.D2")
        
        flog.info('Starting iterative HMM ...')
        heights <- 1:min(min.traverse, 
                         ncol(lesser.frac.data))
        
        boundsnps.pred <- lapply(heights, function(h) {
          
          ct <- cutree(hc, k = h)
          cuts <- unique(ct)
          
          ## look at each group, if deletion present
          boundsnps.pred <- lapply(cuts, function(group) {
            
            if(sum(ct == group)>1) {
              
              mafl <- rowSums(lesser.data[, ct == group]>0)
              sizel <- rowSums(coverage.data[, ct == group]>0)
              
              ## change point
              delta <- c(0, 1)
              z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t),
                                      byrow=TRUE, nrow=2),
                         delta, "binom", list(prob=c(pd, pn)),
                         list(size=sizel), discrete=TRUE)
              results <- Viterbi(z)
              
              ## Get boundaries from states
              boundsnps <- rownames(lesser.frac.data)[results == 1]
              return(boundsnps)
            }
          })
        })
        
        boundsnps_res <- table(unlist(boundsnps.pred))
        
        # vote
        vote <- rep(0, nrow(lesser.frac.data))
        names(vote) <- rownames(lesser.frac.data)
        vote[names(boundsnps_res)] <- boundsnps_res
        
        if(max(vote) == 0) {
          flog.info('Exiting; no new bound SNPs found ...')
          return() ## exit iteration, no more bound SNPs found
        }
        
        vote[vote > 0] <- 1
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.snps.cont <- rep(0, length(vote))
        names(bound.snps.cont) <- names(vote)
        
        for(i in 2:length(vote)) {
          if(vote[i] >= mv & vote[i] == vote[i-1]) {
            bound.snps.cont[i] <- cs
          } else {
            cs <- cs + 1
          }
        }
        
        tb <- table(bound.snps.cont)
        tbv <- as.vector(tb)
        names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0
        
        ## all detected deletions have fewer than 5 SNPs...reached the end
        tbv <- tbv[tbv >= min.num.snps]
        if(length(tbv)==0) {
          flog.info(sprintf('Exiting; less than %s new bound SNPs found ...', min.num.snps))
          return()
        }
        
        HMM_info <- lapply(names(tbv), function(ti) {
          
          bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
          
          ## trim
          bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
          cell_index <<- c(cell_index, list(tumor_subcluster_cells_idx))
          
          return(bound.snps.new)
          
        })
        HMM_output <<- c(HMM_output, HMM_info)
        HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
        
        flog.info(sprintf("Done extracting HMM regions at %s in %s cells ...", 
                          chr, status))
        print(HMM_region)
        
        if(!is.null(HMM_region)){
          hmm.allele.data[infercnv_allele_obj@SNP_info %over% HMM_region,
                          tumor_subcluster_cells_idx] <<- -1
        }
      }
    })
  })
  
  infercnv_allele_obj@expr.data <- hmm.allele.data
  return(infercnv_allele_obj)
  
}
allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples_mod <- function(infercnv_allele_obj,
                                                                      t = 1e-6, pd = 0.1, pn = 0.45,
                                                                      min.num.snps = 5, trim = 0.1,
                                                                      min.traverse = 3){
  
  ## pre-check for allele data
  if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
    flog.info("Initializing the lesser allele fraction ...")
    infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
  }
  
  flog.info("predict_allele_CNV_via_HMM_on_whole_tumor_samples_mod")
  
  chrs = unique(infercnv_allele_obj@gene_order$chr)
  gene_order = infercnv_allele_obj@gene_order
  lesser.frac.data <- infercnv_allele_obj@expr.data
  lesser.data <- infercnv_allele_obj@count.data
  coverage.data <- infercnv_allele_obj@coverage.data
  
  ## initialize hmm allele data
  hmm.allele.data <- matrix(0,
                            nrow = nrow(lesser.data),
                            ncol = ncol(lesser.data))
  rownames(hmm.allele.data) <- rownames(lesser.data)
  colnames(hmm.allele.data) <- colnames(lesser.data)
  
  ## add the normals, so they get predictions too:
  tumor_samples <- c(infercnv_allele_obj@observation_grouped_cell_indices, 
                     infercnv_allele_obj@reference_grouped_cell_indices)
  HMM_output <- c()
  cell_index <- c()
  
  ## run HMM across chromosomes
  lapply(chrs, function(chr){
    
    chr_snp_idx = which(gene_order$chr == chr)
    
    lapply(tumor_samples, function(tumor_sample_cells_idx){
      
      if(length(tumor_sample_cells_idx) > 1){
        
        if(mean(tumor_sample_cells_idx %in% 
                unlist(infercnv_allele_obj@observation_grouped_cell_indices)) == 1){
          status <- "tumor"
        }
        if(mean(tumor_sample_cells_idx %in% 
                unlist(infercnv_allele_obj@reference_grouped_cell_indices)) == 1){
          status <- "normal"
        }
        
        lesser.frac.data <- lesser.frac.data[chr_snp_idx,
                                             tumor_sample_cells_idx,drop=FALSE]
        lesser.data <- lesser.data[chr_snp_idx,
                                   tumor_sample_cells_idx,drop=FALSE]
        coverage.data <- coverage.data[chr_snp_idx,
                                       tumor_sample_cells_idx,drop=FALSE]
        
        d <- parallelDist(t(lesser.frac.data), 
                     method = "euclidean", threads = 5) # parallel dist
        hc <- hclust(d, method="ward.D2")

        flog.info('Starting iterative HMM ...')
        heights <- 1:min(min.traverse, 
                         ncol(lesser.frac.data))

        boundsnps.pred <- lapply(heights, function(h) {
          
          ct <- cutree(hc, k = h)
          cuts <- unique(ct)

          ## look at each group, if deletion present
          boundsnps.pred <- lapply(cuts, function(group) {
            
            if(sum(ct == group)>1) {
              
              mafl <- rowSums(lesser.data[, ct == group]>0)
              sizel <- rowSums(coverage.data[, ct == group]>0)

              ## change point
              delta <- c(0, 1)
              z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t),
                                      byrow=TRUE, nrow=2),
                         delta, "binom", list(prob=c(pd, pn)),
                         list(size=sizel), discrete=TRUE)
              results <- Viterbi(z)

              ## Get boundaries from states
              boundsnps <- rownames(lesser.frac.data)[results == 1]
              return(boundsnps)
            }
          })
        })
        
        boundsnps_res <- table(unlist(boundsnps.pred))
        
        # vote
        vote <- rep(0, nrow(lesser.frac.data))
        names(vote) <- rownames(lesser.frac.data)
        vote[names(boundsnps_res)] <- boundsnps_res

        if(max(vote) == 0) {
          flog.info('Exiting; no new bound SNPs found ...')
          return() ## exit iteration, no more bound SNPs found
        }

        vote[vote > 0] <- 1
        mv <- 1 ## at least 1 vote
        cs <- 1
        bound.snps.cont <- rep(0, length(vote))
        names(bound.snps.cont) <- names(vote)

        for(i in 2:length(vote)) {
          if(vote[i] >= mv & vote[i] == vote[i-1]) {
            bound.snps.cont[i] <- cs
          } else {
            cs <- cs + 1
          }
        }

        tb <- table(bound.snps.cont)
        tbv <- as.vector(tb)
        names(tbv) <- names(tb)
        tbv <- tbv[-1] # get rid of 0

        ## all detected deletions have fewer than 5 SNPs...reached the end
        tbv <- tbv[tbv >= min.num.snps]
        if(length(tbv)==0) {
          flog.info(sprintf('Exiting; less than %s new bound SNPs found ...', min.num.snps))
          return()
        }
        
        HMM_info <- lapply(names(tbv), function(ti) {
          
          bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
          
          ## trim
          bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
          
          cell_index <<- c(cell_index, list(tumor_sample_cells_idx))
          return(bound.snps.new)
        })
        HMM_output <<- c(HMM_output, HMM_info)
        HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
        
        flog.info(sprintf("Done extracting HMM regions at %s in %s cells ...", 
                          chr, status))
        print(HMM_region)
        
        if(!is.null(HMM_region)){
          hmm.allele.data[infercnv_allele_obj@SNP_info %over% HMM_region,
                          tumor_sample_cells_idx] <<- -1
        }
      }
    })
  })
  
  infercnv_allele_obj@expr.data <- hmm.allele.data
  return(infercnv_allele_obj)
  
}
##################


Viterbi.dthmm.allele.adj <- function(object, ...){
  #browser()
  x <- object$x
  size_l <- object$pn$size
  
  if (length(x) < 2) {
    ## not enough run a trace on
    return(2); # neutral state
  }
  
  n <- length(x)
  m <- nrow(object$Pi) # transition matrix
  nu <- matrix(NA, nrow = n, ncol = m)  # scoring matrix
  y <- rep(NA, n) # final trace
  
  ##
  emissions <- matrix(NA, nrow = n, ncol = m) 
  ##
  
  ## init first row
  emission <- dbinom(x = x[1], size = size_l[1], prob = object$pm$prob)
  emission[emission < object$Pi[1,2]] <- object$Pi[1,2]
  ##
  
  ##
  #emission <- 1 / (-1 * emission)
  emission <- emission / sum(emission)
  emissions[1,] <- log(emission)
  #emissions[1,] <- emission
  ##
  
  nu[1, ] <- log(object$delta) + # start probabilities
    emissions[1,]
  
  logPi <- log(object$Pi) # convert transition matrix to log(p)
  
  for (i in 2:n) {
    
    matrixnu <- matrix(nu[i - 1, ], nrow = m, ncol = m)
    
    ##
    emission <- dbinom(x = x[i], size = size_l[i], prob = object$pm$prob)
    emission[emission < object$Pi[1,2]] <- object$Pi[1,2]
    
    #emission <- 1 / (-1 * emission)
    emission <- emission / sum(emission)
    emissions[i,] <- log(emission)
    #emissions[i,] <- emission
    ##
    
    nu[i, ] <- apply(matrixnu + logPi, 2, max) + emissions[i, ] 
    
  }
  #return(nu)
  
  if (any(nu[n, ] == -Inf)) 
    stop("Problems With Underflow")
  
  ## traceback
  y[n] <- which.max(nu[n, ])
  
  for (i in seq(n - 1, 1, -1))
    y[i] <- which.max(logPi[, y[i + 1]] + nu[i, ])
  
  return(y)
}

## aggregrate into a higher level
# agg_method <- function(data, index, method = c("mean","median","sum")){
#   
#   method <- match.arg(method)
#   
#   tmp <- lapply(seq_len(ncol(data)), 
#                 function(x){
#                   tapply(data[,x], 
#                          index, 
#                          method)
#                 })
#   tmp <- do.call(cbind, tmp)
#   colnames(tmp) <- colnames(data)
#   tmp <- round(tmp)
#   
#   return(tmp)
#   
# }

## aggregate HMM-based obj (snp level) into HMM_based obj (gene level)
# aggregate_gene <- function(infercnv_allele_obj, gene_annot,
#                            state_agg_method = c("mean","median"),
#                            snp_agg_method = c("mean","median", "sum")){
#   state_agg_method <- match.arg(state_agg_method)
#   snp_agg_method <- match.arg(snp_agg_method)
#   
#   flog.info(sprintf("Using %s method to aggregrate snp-based HMM into gene-based HMM", 
#                     state_agg_method))
#   gene_state <- agg_method(infercnv_allele_obj@expr.data,
#                            infercnv_allele_obj@SNP_info$gene,
#                            state_agg_method)
#   
#   flog.info(sprintf("Using %s method to aggregrate snp into gene level", 
#                     snp_agg_method))
#   gene_count.data <- agg_method(infercnv_allele_obj@count.data,
#                                 infercnv_allele_obj@SNP_info$gene,
#                                 snp_agg_method)
#   gene_allele.data <- agg_method(infercnv_allele_obj@allele.data,
#                                  infercnv_allele_obj@SNP_info$gene,
#                                  snp_agg_method)
#   gene_coverage.data <- agg_method(infercnv_allele_obj@coverage.data,
#                                    infercnv_allele_obj@SNP_info$gene,
#                                    snp_agg_method)
#   
#   gene_order <- gene_annot[rownames(gene_annot) %in% rownames(gene_state),]
#   gene_state <- gene_state[rownames(gene_order), ]
#   gene_count.data <- gene_count.data[rownames(gene_order), ]
#   gene_allele.data <- gene_allele.data[rownames(gene_order), ]
#   gene_coverage.data <- gene_coverage.data[rownames(gene_order), ]
#   SNP_Gene <- GenomicRanges::GRanges(gene_order[[C_CHR]],
#                                      IRanges::IRanges(gene_order[[C_START]],
#                                                       gene_order[[C_STOP]]))
#   names(SNP_Gene) <- rownames(gene_order)
#   
#   infercnv_allele_obj@expr.data <- gene_state
#   infercnv_allele_obj@count.data <- gene_count.data
#   infercnv_allele_obj@allele.data <- gene_allele.data
#   infercnv_allele_obj@coverage.data <- gene_coverage.data
#   infercnv_allele_obj@SNP_info <- SNP_Gene
#   infercnv_allele_obj@gene_order <- gene_order
#   
#   validate_infercnv_allele_obj(infercnv_allele_obj)
#   
#   return(infercnv_allele_obj)
# }

## @title calAlleleBoundaries
## 
## @description Determines the LOHs/Deletions boundaries using HMM model
## 
## @param infercnv_allele_obj infercnv_allele object
## 
## @param distance_method A method used for pre-processing NA value
## 
## @param ncores Number of cores used for parallel task
## 
## @param min.traverse The max number of groups used for cutting tree, default 3
## 
## @param t Transition probability, default 1e-6
## 
## @param pd Emission probability for LOH/deletion
## 
## @param Pn Emission probability for neutral
## 
## @param min.num.snps A cutoff of minimal length of SNPs used for determining LOH/deletion
## 
## @param trim
## 
## @return A list containing putative LOH/deletion boundary
## 
## @importFrom parallelDist parDist
## 
## @keywords internal
## 
## @noRd
# calAlleleBoundaries <- function(infercnv_allele_obj, mode = c("all", "tumor", "normal"),
#                                 distance_method = c("Filter_threshold", "Remove_NA", "HB"), ncores = 20,
#                                 min.traverse = 3, t = 1e-6, pd = 0.1, pn = 0.45, 
#                                 min.num.snps = 5, trim = 0.1){
#   
#   if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
#     flog.info("Initializing the lesser allele fraction ...")
#     infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
#   }
#   
#   ###
#   analysis_mode <- match.arg(mode)
#   flog.info(sprintf("Using %s mode ...", analysis_mode))
#   
#   fraction_method <- match.arg(distance_method)
#   flog.info(sprintf("Using %s to process matrix ...", fraction_method))
#   ###
#   allele.lesser.data <- infercnv_allele_obj@count.data
#   allele.data <- infercnv_allele_obj@allele.data
#   allele.coverage.data <- infercnv_allele_obj@coverage.data
#   
#   ###
#   if (analysis_mode == "tumor"){
#     allele.lesser.data <- allele.lesser.data[, unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
#     allele.data <- allele.data[, unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
#     allele.coverage.data <- allele.coverage.data[, unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
#     
#     mat.tot <- infercnv_allele_obj@expr.data[, unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
#   } else if (analysis_mode == "normal"){
#     allele.lesser.data <- allele.lesser.data[, unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
#     allele.data <- allele.data[, unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
#     allele.coverage.data <- allele.coverage.data[, unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
#     
#     mat.tot <- infercnv_allele_obj@expr.data[, unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
#   } else{
#     mat.tot <- infercnv_allele_obj@expr.data
#   }
#   ###
#   flog.info("Starting cluster cells from population ...")
#   
#   if(fraction_method == "Filter_threshold"){
#     mat.tot[is.na(mat.tot)] <- 0 # omit no coverage
#     mat.tot[allele.data == 0 & allele.coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
#     mat.tot[allele.coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
#     
#     mat.smooth <- apply(mat.tot, 2, runmean, k=31)
#     
#     flog.info("Starting calculate distance ...")
#     d <- parDist(t(mat.smooth), method = "euclidean", threads = ncores) # parallel dist
#   }
#   else if(fraction_method == "Remove_NA"){
#     mat.smooth <- apply(mat.tot, 2, runmean, k=31)
#     mat.smooth[is.na(mat.smooth)] <- 0
#     
#     flog.info("Starting calculate distance ...")
#     d <- parDist(t(mat.smooth), method = "euclidean", threads = ncores) # parallel dist
#   }
#   else if(fraction_method == "HB"){
#     mat.smooth <- apply(mat.tot, 2, runmean, k=31)
#     flog.info("Starting calculate distance ...")
#     d <- dist(t(mat.smooth), method = "euclidean") # too slow
#     d[is.na(d)] <- 0
#     d[is.nan(d)] <- 0
#     d[is.infinite(d)] <- 0
#   }
#   
#   flog.info("Starting calculate Hierarchical Clustering ...")
#   hc <- hclust(d, method="ward.D2")
#   
#   flog.info('Starting iterative HMM ...')
#   heights <- 1:min(min.traverse, ncol(allele.lesser.data))
#   
#   boundsnps.pred <- lapply(heights, function(h) {
#     
#     ct <- cutree(hc, k = h)
#     cuts <- unique(ct)
#     
#     ## look at each group, if deletion present
#     boundsnps.pred <- lapply(cuts, function(group) {
#       
#       if(sum(ct == group)>1) {
#         mafl <- rowSums(allele.lesser.data[, ct == group]>0)
#         sizel <- rowSums(allele.coverage.data[, ct == group]>0)
#         
#         ## change point
#         delta <- c(0, 1)
#         z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
#                                 byrow=TRUE, nrow=2), 
#                    delta, "binom", list(prob=c(pd, pn)), 
#                    list(size=sizel), discrete=TRUE)
#         results <- Viterbi(z)
#         
#         ## Get boundaries from states
#         boundsnps <- rownames(allele.lesser.data)[results == 1]
#         return(boundsnps)
#         
#       }
#     })
#   })
#   
#   boundsnps_res <- table(unlist(boundsnps.pred))
#   
#   ## vote
#   vote <- rep(0, nrow(allele.lesser.data))
#   names(vote) <- rownames(allele.lesser.data)
#   vote[names(boundsnps_res)] <- boundsnps_res
#   
#   if(max(vote) == 0) {
#     flog.info('Exiting; no new bound SNPs found ...')
#     return() ## exit iteration, no more bound SNPs found
#   }
#   
#   vote[vote > 0] <- 1
#   mv <- 1 ## at least 1 vote
#   cs <- 1
#   bound.snps.cont <- rep(0, length(vote))
#   names(bound.snps.cont) <- names(vote)
#   
#   for(i in 2:length(vote)) {
#     if(vote[i] >= mv & vote[i] == vote[i-1]) {
#       bound.snps.cont[i] <- cs
#     } else {
#       cs <- cs + 1
#     }
#   }
#   
#   tb <- table(bound.snps.cont)
#   tbv <- as.vector(tb)
#   names(tbv) <- names(tb)
#   tbv <- tbv[-1] # get rid of 0
#   
#   ## all detected deletions have fewer than 5 SNPs...reached the end
#   tbv <- tbv[tbv >= min.num.snps]
#   if(length(tbv)==0) {
#     flog.info(sprintf('Exiting; less than %s new bound SNPs found ...', min.num.snps))
#     return()
#   }
#   
#   HMM_info <- lapply(names(tbv), function(ti) {
#     
#     bound.snps.new <- names(bound.snps.cont)[bound.snps.cont == ti]
#     
#     ## trim
#     bound.snps.new <- bound.snps.new[1:round(length(bound.snps.new)-length(bound.snps.new)*trim)]
#     
#     return(bound.snps.new)
#     
#   })
#   HMM_region <- do.call("c", lapply(HMM_info, function(bs) range(infercnv_allele_obj@SNP_info[bs])))
#   
#   flog.info("Done extracting HMM regions ...")
#   
#   return(list("HMM_info" = HMM_info,
#               "HMM_region" = HMM_region))
# }