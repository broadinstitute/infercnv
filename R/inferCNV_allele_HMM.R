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
        
        mafl <- rowSums(lesser.data[chr_snp_idx, tumor_subcluster_cells_idx,drop=FALSE] > 0)
        sizel <- rowSums(coverage.data[chr_snp_idx, tumor_subcluster_cells_idx,drop=FALSE] > 0)
        
        ## change point
        delta <- c(0, 1)
        z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                byrow=TRUE, nrow=2), 
                   delta, "binom", list(prob=c(pd, pn)), 
                   list(size=sizel), discrete=TRUE)
        results <- Viterbi(z)
        
        boundsnps <- rownames(lesser.data[chr_snp_idx,tumor_subcluster_cells_idx,drop=FALSE])[results == 1]
        
        ## vote
        vote <- rep(0, nrow(lesser.data[chr_snp_idx,tumor_subcluster_cells_idx,drop=FALSE]))
        names(vote) <- rownames(lesser.data[chr_snp_idx,tumor_subcluster_cells_idx,drop=FALSE])
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
          
          return(bound.snps.new)
          
        })
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
        
        mafl <- rowSums(lesser.data[chr_snp_idx, tumor_sample_cells_idx,drop=FALSE] > 0)
        sizel <- rowSums(coverage.data[chr_snp_idx, tumor_sample_cells_idx,drop=FALSE] > 0)
        
        ## change point
        delta <- c(0, 1)
        z <- dthmm(mafl, matrix(c(1-t, t, t, 1-t), 
                                byrow=TRUE, nrow=2), 
                   delta, "binom", list(prob=c(pd, pn)), 
                   list(size=sizel), discrete=TRUE)
        results <- Viterbi(z)
        
        boundsnps <- rownames(lesser.data[chr_snp_idx,tumor_sample_cells_idx,drop=FALSE])[results == 1]
        
        ## vote
        vote <- rep(0, nrow(lesser.data[chr_snp_idx,tumor_sample_cells_idx,drop=FALSE]))
        names(vote) <- rownames(lesser.data[chr_snp_idx,tumor_sample_cells_idx,drop=FALSE])
        
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
          
          return(bound.snps.new)
          
        })
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

#' @title plot_allele
#' 
#' @description plot a summary figure containing allele frequency, HMM prediction 
#' 
#' @keywords internal
#' 
#' @noRd
plot_allele <- function(infercnv_allele_obj, 
                        #expression_mode = F,
                        #allele_frequency_mode = F,
                        #HMM = NULL,
                        #use_common_gene = F, 
                        name_to_plot,
                        trend_smK = 31,
                        CELL_POINT_ALPHA = 0.6, dotsize=0.3, colorscheme = "BlueRed"){
  
  options(bitmapType = "cairo")
  ############# cell annotation
  normal_cells = colnames(infercnv_allele_obj@allele.data)[unlist(infercnv_allele_obj@reference_grouped_cell_indices)]
  malignant_cells = colnames(infercnv_allele_obj@allele.data)[unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
  
  num_normal_cells = length(normal_cells)
  #############
  
  ############# allele data loading
  # if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
  #   flog.info("Initializing the lesser allele fraction ...")
  #   infercnv_allele_obj <- infercnv:::setAlleleMatrix(infercnv_allele_obj)
  # }
  # allele_matrix <- infercnv_allele_obj@expr.data
  # #############
  
  ############# allele data
  flog.info("Building snp plot ...")
  
  allele_matrix <- infercnv_allele_obj@allele.data/infercnv_allele_obj@coverage.data # fraction
  #allele_matrix <- infercnv_allele_obj@allele.lesser.data/infercnv_allele_obj@coverage.data

  allele_matrix[is.na(allele_matrix)] <- 0 # omit no coverage
  allele_matrix[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
  allele_matrix[infercnv_allele_obj@coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
  
  min.cells = 3
  cells_w_ref_allele = rowSums(allele_matrix != 0 & allele_matrix < 0.5)
  cells_w_alt_allele = rowSums(allele_matrix != 0 & allele_matrix > 0.5)
  
  allele_matrix = allele_matrix[(cells_w_ref_allele >= min.cells & cells_w_alt_allele >= min.cells), ]
  num_snps_all = nrow(allele_matrix)
  flog.info(sprintf("Number of heterozygous snps used for plotting: %s ...", num_snps_all))
  
  flog.info("Setting alt allele fraction to the tumor cell-population minor allele ...")
  
  mAF_allele_matrix = apply(allele_matrix, 1, function(x) {
    nonzero_val_idx = which(x>0)
    nonzero_vals = x[nonzero_val_idx]
    
    frac_high = sum(nonzero_vals>0.5)/length(nonzero_vals)
    
    ## focus allele selection based on the tumor cells only.
    tumor_vals = x[unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
    tumor_nonzero_vals = tumor_vals[tumor_vals>0]
    if (length(tumor_nonzero_vals) > 0) {
      frac_high = sum(tumor_nonzero_vals>0.5)/length(tumor_nonzero_vals)
    }
    
    if ( frac_high > 0.5) {
      x[x==1] = 0.999
      x[nonzero_val_idx ] = 1 - x[nonzero_val_idx]
    }
    x
  })
  
  allele_matrix = t(mAF_allele_matrix)
  
  if (!is.null(infercnv_allele_obj@tumor_subclusters)) {
    ## define cell ordering.
    flog.info("Exracting clustering info ...")
    ordered_cells <- colnames(allele_matrix)[unlist(infercnv_allele_obj@tumor_subclusters$subclusters)]
  }
  
  flog.info("Melting matrix ...")
  alleledatamelt = melt(as.matrix(allele_matrix))
  colnames(alleledatamelt) = c('chrpos', 'cell', 'AF')
  
  if (!is.null(infercnv_allele_obj@tumor_subclusters)) {
    alleledatamelt$cell = ordered(alleledatamelt$cell, levels=ordered_cells)
  }
  
  alleledatamelt = alleledatamelt %>% separate(chrpos, ":", into=c('seqnames', 'pos', "end"), remove=FALSE)
  
  alleledatamelt$chr = str_replace(string=alleledatamelt$seqnames, pattern="chr", replacement="")
  
  alleledatamelt = alleledatamelt %>% dplyr::filter(chr %in% 1:22)
  
  alleledatamelt = alleledatamelt %>% mutate(chr = ordered(chr, levels=1:22))
  alleledatamelt$pos = as.numeric(alleledatamelt$pos)
  
  ## get chr bounds for plotting later.
  chr_maxpos_snp = alleledatamelt %>% group_by(chr) %>% summarise(maxpos = max(pos))
  chr_maxpos_snp$minpos = 1
  
  
  alleledatamelt = alleledatamelt %>% dplyr::filter(AF > 0)  ## if AF==0, then means we have no coverage.
  
  alleledatamelt = alleledatamelt %>% mutate(cellchr = paste(cell, seqnames, sep=":"))
  
  midpt = mean(alleledatamelt$AF)
  
  alleledatamelt$sample_type = "tumor"
  if (num_normal_cells > 0) {
    alleledatamelt$sample_type[ alleledatamelt$cell %in% normal_cells ] = "normal"
    alleledatamelt$sample_type <- factor(alleledatamelt$sample_type)
  }
  #############
  
  ############# expression data
  gexdata <- NULL
  chr_maxpos_gex <- NULL
  # if(expression_mode){
  #   
  #   if(num_normal_cells == 0){
  #     flog.info("plot_allele:: No reference here")
  #   }
  #   
  #   ## annotation preparation
  #   gencode_gene_pos <- infercnv_allele_obj@gene_order
  #   gencode_gene_pos$chr = str_replace(string=gencode_gene_pos$chr, pattern="chr", replacement="")
  #   gencode_gene_pos$genename <- gencode_gene_pos %>% row.names()
  #   #gencode_gene_pos = gencode_gene_pos %>% dplyr::filter(chr %in% 1:22 & genename %in% infercnv_allele_obj@SNP_info$gene_name) ## modification
  #   gencode_gene_pos = gencode_gene_pos %>% dplyr::filter(chr %in% 1:22)
  #   gencode_gene_pos = gencode_gene_pos %>% mutate(chr = ordered(chr, levels=1:22))
  #   gencode_gene_pos$pos = as.numeric( (gencode_gene_pos$start + gencode_gene_pos$stop)/2)
  #   
  #   chr_maxpos_gex = gencode_gene_pos %>% group_by(chr) %>% summarise(maxpos = max(pos))
  #   chr_maxpos_gex$minpos = 1
  #   ##
  #   #gex_malignant_cells <- colnames(express_matrix)[colnames(express_matrix) %in% malignant_cells]
  #   #gex_normal_cells <- colnames(express_matrix)[colnames(express_matrix) %in% normal_cells]
  #   #express_matrix <- express_matrix[,c(gex_malignant_cells, gex_normal_cells)]
  #   express_matrix <- infercnv_allele_obj@expr.data
  #   if(num_normal_cells > 0){
  #     express_index <- (rowMeans(express_matrix[, malignant_cells]) > 1 & 
  #                         rowMeans(express_matrix[, normal_cells]) > 1)
  #     cell_index <- apply(express_matrix[, malignant_cells], 1, function(x) { sum(x>0 & ! is.na(x)) >= 3}) &
  #       apply(express_matrix[, normal_cells], 1, function(x) { sum(x>0 & ! is.na(x)) >= 3})
  #   }
  #   else{
  #     express_index <- (rowMeans(express_matrix) > 1)
  #     cell_index <- apply(express_matrix, 1, function(x) { sum(x>0 & ! is.na(x)) >= 3})
  #   }
  #   #gene_index <- rownames(express_matrix) %in% infercnv_allele_obj@SNP_info$gene_name
  #   
  #   express_matrix <- express_matrix[express_index & cell_index,]
  #   
  #   gexdata = t( t(express_matrix)/colSums(express_matrix)) * 1e6 # cpm normalization
  #   log2data = log2(gexdata+1)
  #   zscaled_data = t(scale(t(gexdata), scale=T, center=T))
  #   zscaled_data[is.nan(zscaled_data)] = 0
  #   
  #   ## subtract normals:
  #   if(num_normal_cells > 0){
  #     normal_reference_data = zscaled_data[,normal_cells]
  #     gexp.norm = zscaled_data[,malignant_cells]
  #   }
  #   else{
  #     normal_reference_data <- zscaled_data
  #     gexp.norm <- zscaled_data
  #   }
  #   
  #   ## subtract tumor from normal:
  #   normal_gene_means = rowMeans(normal_reference_data)
  #   
  #   for (i in 1:length(normal_gene_means) ) {
  #     gexp.norm[i,] = gexp.norm[i,]  - normal_gene_means[i]
  #   }
  #   
  #   gexp.melt = melt(gexp.norm)
  #   colnames(gexp.melt) = c('genename', 'cell', 'exp.norm')
  #   
  #   gexdata = left_join(gexp.melt, gencode_gene_pos, key='genename')
  #   
  #   gexdata = gexdata %>% dplyr::filter(! is.na(chr))
  #   
  #   ## do smoothing by cell and chr.
  #   
  #   gexdata = gexdata %>% mutate(cellchr = paste(cell, chr, sep=":"))
  #   splitdata = split(gexdata, gexdata$cellchr)
  #   
  #   smoother = function(df) {
  #     df = df %>% arrange(pos)
  #     df$exp.norm.smoothed = runmean(df$exp.norm,k=101, align="center")
  #     ## recenter cells
  #     
  #     return(df)
  #   }
  #   gexdata = do.call(rbind, lapply(splitdata, smoother))
  #   
  #   ## mean center cells again after smoothing
  #   splitdata = split(gexdata, gexdata$cell)
  #   gexdata = do.call(rbind, lapply(splitdata, function(x) {
  #     m = mean(x$exp.norm.smoothed)
  #     x$exp.norm.smoothed = x$exp.norm.smoothed - m
  #     return(x)
  #   } ) )
  #   
  #   ## cluster cells according to expr pattern:
  #   gene_expr_matrix = gexdata %>% select(genename, cell, exp.norm.smoothed) %>% 
  #     spread(key=cell, value=exp.norm.smoothed)
  #   
  #   rownames(gene_expr_matrix) = gene_expr_matrix$genename
  #   gene_expr_matrix = gene_expr_matrix[,-1]
  #   h = hclust(dist(t(gene_expr_matrix)))
  #   ordered_cells = colnames(gene_expr_matrix)[h$order]
  #   
  #   gexdata$cell = ordered(gexdata$cell, levels=ordered_cells)
  #   
  #   ## set up base plot
  #   ## define chr boundaries based on max coordinates for now.
  #   
  #   q = quantile(gexdata$exp.norm.smoothed, c(0.025, 0.975))
  #   
  #   gexdata$exp.norm.smoothed[gexdata$exp.norm.smoothed < q[1] ] = q[1]
  #   gexdata$exp.norm.smoothed[gexdata$exp.norm.smoothed > q[2] ] = q[2]
  #   
  #   noise_range = 0.2
  #   gexdata$exp.norm.smoothed[gexdata$exp.norm.smoothed > -noise_range & gexdata$exp.norm.smoothed < noise_range] = 0
  # }
  # #############
  # 
  
  make_plots = function(alleledataToPlot, gexdataToPlot = NULL, 
                        chr_maxpos_snp = chr_maxpos_snp, chr_maxpos_gex = chr_maxpos_gex,
                        infercnv_allele_obj) {
    #browser()
    ## normal cell plot
    
    normal_snps_plot = NA
    bar_plot = NA
    expr_plot = NA
    malignant_HMM_plot = NA
    color_cell <- "#F8766D"
    
    if (num_normal_cells > 0) {
      
      color_cell <- c("#00BFC4", "#F8766D")
      flog.info("Making normal plot ...")
      
      normal_dataToPlot = alleledataToPlot %>% dplyr::filter(cell %in% normal_cells)
      
      
      normal_snps_plot = ggplot(data=normal_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
        
        labs(title=paste0("SNPs in normal cells: ", length(normal_cells))) +
        theme_bw() +
        theme(text = element_text(size = 12),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank()
        ) +
        
        geom_vline(data=chr_maxpos_snp, aes(xintercept=minpos), color=NA) +
        geom_vline(data=chr_maxpos_snp, aes(xintercept=maxpos), color=NA) +
        
        geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius() 
      
      
      if (colorscheme == "BlueRed") {
        normal_snps_plot = normal_snps_plot + scale_color_gradient(low="blue", high="red")
      } else {
        normal_snps_plot = normal_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
      }
    }
    
    ## malignant cell plot
    
    flog.info("Making malignant plot ...")
    
    num_malignant_cells = length(malignant_cells)
    
    malignant_dataToPlot = alleledataToPlot %>% dplyr::filter(cell %in% malignant_cells)
    
    malignant_snps_plot = ggplot(data=malignant_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
      
      labs(title=paste0("SNPs in tumor cells: ", length(malignant_cells))) +
      theme_bw() +
      theme(text = element_text(size = 12),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      ) +
      
      geom_vline(data=chr_maxpos_snp, aes(xintercept=minpos), color=NA) +
      geom_vline(data=chr_maxpos_snp, aes(xintercept=maxpos), color=NA) +
      
      geom_point(aes(x=pos, y=cell, color=AF), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()
    
    if (colorscheme == "BlueRed") {
      malignant_snps_plot = malignant_snps_plot + scale_color_gradient(low="blue", high="red")
    } else {
      malignant_snps_plot = malignant_snps_plot + scale_color_gradient2(low="blue", mid='yellow', high="red", midpoint=midpt)
    }
    
    
    flog.info("Making allele freq plot ...")
    
    
    allele_freq_means = alleledataToPlot %>%
      group_by(chrpos,sample_type) %>%
      mutate(grp_pos_mean_AF = mean(AF)) %>% select(chrpos, chr, pos, sample_type, grp_pos_mean_AF) %>% unique()
    
    
    
    ## smooth mean AFs across chromosomes for each sample type.
    flog.info("Smoothing sample means for trend lines ...")
    
    allele_freq_means = allele_freq_means %>% mutate(sampleTypeChr = paste(sample_type, chr, sep=":"))
    splitdata = split(allele_freq_means, allele_freq_means$sampleTypeChr)
    
    smoother = function(df) {
      df = df %>% arrange(pos)
      df$grp_pos_mean_AF_sm = runmean(df$grp_pos_mean_AF, k=trend_smK, align="center")
      return(df)
    }
    
    allele_freq_means = do.call(rbind, lapply(splitdata, smoother))
    
    allele_freq_plot = allele_freq_means %>%
      ggplot(aes(x=pos, y=grp_pos_mean_AF)) +
      facet_grid (~chr, scales = 'free_x', space = 'fixed') +
      
      labs(title="Allele frequency",
           y="mean_AF") +
      theme_bw() +
      theme(text = element_text(size = 12),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank()
      ) +
      geom_vline(data=chr_maxpos_snp, aes(xintercept=minpos), color=NA) +
      geom_vline(data=chr_maxpos_snp, aes(xintercept=maxpos), color=NA) +
      geom_point(aes(color=sample_type), alpha=0.2, size=0.2)
    
    allele_freq_plot_w_trendlines = allele_freq_plot +
      geom_line(data=allele_freq_means,
                aes(x=pos, y=grp_pos_mean_AF_sm, color=sample_type), size=0.5, alpha=1) +
      scale_color_manual(values = color_cell)
    
    # # HMM prediction if existed
    # if(!is.null(HMM)){
    #   
    #   flog.info("Making HMM prediction plot ...")
    #   
    #   # malignant_dataToPlot_Gr <- with(malignant_dataToPlot, 
    #   #                                 GenomicRanges::GRanges(seqnames, 
    #   #                                                        IRanges::IRanges(as.numeric(as.character(pos)), 
    #   #                                                                         as.numeric(as.character(end)))))
    #   # malignant_dataToPlot$HMM_region <- "Neural"
    #   # malignant_dataToPlot$HMM_region[malignant_dataToPlot_Gr %over% infercnv_allele_obj@HMM_region] <- "Candidate Deletions/LOHs"
    #   # 
    #   # malignant_HMM_plot = ggplot(data=malignant_dataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
    #   #   
    #   #   labs(title="HMM prediction plot") +
    #   #   theme_bw() +
    #   #   theme(text = element_text(size = 6),
    #   #         axis.ticks.x = element_blank(),
    #   #         axis.text.x = element_blank(),
    #   #         axis.title.x = element_blank(),
    #   #         axis.ticks.y = element_blank(),
    #   #         axis.text.y = element_blank(),
    #   #         panel.grid.major.x = element_blank(),
    #   #         panel.grid.minor.x = element_blank(),
    #   #         panel.grid.major.y = element_blank(),
    #   #         panel.grid.minor.y = element_blank()
    #   #   ) +
    #   #   geom_vline(data=chr_maxpos_snp, aes(xintercept=minpos), color=NA) +
    #   #   geom_vline(data=chr_maxpos_snp, aes(xintercept=maxpos), color=NA) +
    #   #   geom_point(aes(x=pos, y=cell, color=HMM_region), alpha=CELL_POINT_ALPHA, size=dotsize) + scale_radius()
    #   tmp_data <- HMM %>% as.data.frame()
    #   tmp_data <- tmp_data[tmp_data$seqnames %in% paste0("chr",1:22),]
    #   tmp_data$chr <- str_replace(string=tmp_data$seqnames, pattern="chr", replacement="") %>% 
    #     ordered(levels = 1:22)
    #   
    #   malignant_HMM_plot <- ggplot(data=tmp_data) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
    #     
    #     labs(title="LOHs/Deletions boundary (HMM model)") +
    #     theme_bw() +
    #     theme(text = element_text(size = 12),
    #           axis.ticks.x = element_blank(),
    #           axis.text.x = element_blank(),
    #           axis.title.x = element_blank(),
    #           axis.ticks.y = element_blank(),
    #           axis.text.y = element_blank(),
    #           panel.grid.major.x = element_blank(),
    #           panel.grid.minor.x = element_blank(),
    #           panel.grid.major.y = element_blank(),
    #           panel.grid.minor.y = element_blank()
    #     ) +
    #     geom_vline(data=chr_maxpos_snp, aes(xintercept=minpos), color=NA) +
    #     geom_vline(data=chr_maxpos_snp, aes(xintercept=maxpos), color=NA) +
    #     geom_rect(data = tmp_data, aes(xmin = start, xmax = end, ymin=0, ymax =1,
    #                                    alpha = 0.8), 
    #               show.legend = T, fill = "blue")
    # }
    # 
    # if(!is.null(gexdataToPlot)){
    #   
    #   flog.info("Making expression plot")
    #   expr_plot = ggplot(data=gexdataToPlot) + facet_grid (~chr, scales = 'free_x', space = 'fixed') +
    #     
    #     labs(title="Gene expression") +
    #     theme_bw() +
    #     theme(text = element_text(size = 12),
    #           axis.ticks.x = element_blank(),
    #           axis.text.x = element_blank(),
    #           axis.title.x = element_blank(),
    #           axis.ticks.y = element_blank(),
    #           axis.text.y = element_blank(),
    #           panel.grid.major.x = element_blank(),
    #           panel.grid.minor.x = element_blank(),
    #           panel.grid.major.y = element_blank(),
    #           panel.grid.minor.y = element_blank()
    #     ) +
    #     
    #     geom_vline(data=chr_maxpos_gex, aes(xintercept=minpos), color=NA) +
    #     geom_vline(data=chr_maxpos_gex, aes(xintercept=maxpos), color=NA) +
    #     
    #     geom_point(aes(x=pos, y=cell, color=exp.norm.smoothed), alpha=0.6, size=dotsize) +
    #     
    #     scale_colour_gradient2(low = "blue",
    #                            mid = "white",
    #                            high = "red",
    #                            midpoint = 0,
    #                            space = "Lab",
    #                            guide = "colourbar",
    #                            aesthetics = "colour")
    # }
    
    # if(allele_frequency_mode){
    #   
    #   cell_count <- alleledatamelt %>% select(chr, cell, sample_type) %>% unique() %>% 
    #     group_by(chr, sample_type) %>% summarise(cell_count = n(), .groups = 'drop') %>% 
    #     select(cell_count)
    #   bar_data <- alleledatamelt %>% select(chrpos, chr, sample_type) %>% unique() %>% 
    #     group_by(chr, sample_type) %>% summarise(snp_count = n(), .groups = 'drop') %>% 
    #     mutate(cell_count, "SNPs coverage per cell" = snp_count/cell_count)
    #   
    #   bar_plot <- bar_data %>% ggplot(aes(x = sample_type, y = `SNPs coverage per cell`, fill = sample_type)) +
    #     geom_bar(stat='identity') +
    #     facet_grid (~chr, scales = 'free_x', space = 'fixed') +
    #     
    #     labs(title="Fraction of heterozygous snps used") +
    #     theme_bw() +
    #     theme(text = element_text(size = 12),
    #           axis.ticks.x = element_blank(),
    #           axis.text.x = element_blank(),
    #           axis.title.x = element_blank(),
    #           panel.grid.major.x = element_blank(),
    #           panel.grid.minor.x = element_blank(),
    #           panel.grid.major.y = element_blank(),
    #           panel.grid.minor.y = element_blank()
    #     ) +
    #     scale_fill_manual(values = color_cell)
    #   
    # }
    #     
    #ratio_normal_cells = max(0.25, num_normal_cells/(num_normal_cells + num_malignant_cells))
    # if(is.null(expr_plot)){
    #   if(num_normal_cells > 0) {
    #     pg = plot_grid(normal_snps_plot, allele_freq_plot_w_trendlines, bar_plot,
    #                    malignant_HMM_plot, malignant_snps_plot,
    #                    ncol=1, align='v', rel_heights=rep(1/5,5))
    #     return(pg)
    #   } else{
    #     pg = plot_grid(allele_freq_plot_w_trendlines, bar_plot,
    #                    malignant_HMM_plot, malignant_snps_plot, 
    #                    ncol=1, align='v', rel_heights=rep(1/4,4))
    #     return(pg)
    #   }
    #   
    # } else{
    #   if(num_normal_cells > 0){
    #     pg = plot_grid(normal_snps_plot, allele_freq_plot_w_trendlines,
    #                    malignant_HMM_plot, malignant_snps_plot, expr_plot,
    #                    ncol=1, align='v', rel_heights=rep(1/5,5))
    #     return(pg)
    #   } else{
    #     pg = plot_grid(allele_freq_plot_w_trendlines, bar_plot,
    #                    malignant_HMM_plot, malignant_snps_plot, expr_plot,
    #                    ncol=1, align='v', rel_heights=rep(1/5,5))
    #     return(pg)
    #   }
    # }
    
    if (num_normal_cells > 0) {
      ratio_normal_cells = max(0.25, num_normal_cells/(num_normal_cells + num_malignant_cells))
      pg = plot_grid(normal_snps_plot, 
                     allele_freq_plot_w_trendlines, 
                     malignant_snps_plot, 
                     ncol=1, align='v', 
                     rel_heights=c(ratio_normal_cells, 0.25, 1-ratio_normal_cells))
    } else{
      pg = plot_grid(allele_freq_plot_w_trendlines, 
                     malignant_snps_plot,
                     ncol=1, align='v',
                     rel_heights=c(1/3,2/3))
    }
    return(pg)
  }
  
  flog.info("Generating outputs ...")

  # if(use_common_gene & expression_mode){
  #   tmp_snp <- infercnv_allele_obj@SNP_info[unique(alleledatamelt$chrpos)]
  #   common_gene <- intersect(tmp_snp$gene_name,
  #                            unique(gexdata$genename))
  #   snp_index <- tmp_snp[tmp_snp$gene_name %in% common_gene] %>% names()
  #   
  #   alleledatamelt <- alleledatamelt %>% dplyr::filter(chrpos %in% snp_index)
  #   gexdata <- gexdata %>% dplyr::filter(genename %in% common_gene)
  # }
  # 
  p <- make_plots(alleledataToPlot = alleledatamelt, 
                  gexdataToPlot = gexdata,
                  chr_maxpos_snp = chr_maxpos_snp, chr_maxpos_gex = chr_maxpos_gex, 
                  infercnv_allele_obj = infercnv_allele_obj)
  
  ggsave (name_to_plot, p, width = 13.33, height = 7.5, units = 'in', dpi = 300)
  flog.info("Done!")
}
