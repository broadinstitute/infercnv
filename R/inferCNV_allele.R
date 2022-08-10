#' The infercnv_allele class
#' 
#' @description This class extends the functionality of infercnv class aiming to
#' incorporate the allele information.
#' 
#' Slots in the infercnv_allele object (besides infercnv) include:
#' 
#' @slot expr.data <matrix> the lesser allele fraction matrix (infercnv_allele obj).
#' 
#' @slot count.data <matrix> the lesser allele count matrix (infercnv_allele obj).
#'
#' @slot allele.data <matrix> the alternate allele data matrix.
#' 
#' @slot coverage.data <matrix> the coverage allele data matrix.
#' 
#' @slot SNP_info <GRanges> A GRanges object stores the SNP sites information given allele data.
#' 
#' @export
infercnv_allele <- methods::setClass(
  "infercnv_allele",
  slots = c(
    expr.data = "ANY",
    count.data = "ANY",
    allele.data = "ANY",
    coverage.data = "ANY",
    SNP_info = "GRanges"),
  contains = "infercnv")  

#' @title validate_infercnv_allele_obj()
#'
#' @description validate an infercnv_allele_obj
#' ensures that allele data and coverage data do have the same cell population,
#' and indicate the same SNP region as @SNP_info stored.
#' Otherwise, throws an error and stops execution.
#'
#' @param infercnv_allele_obj infercnv_allele_object
#'
#' @return none
#' 
#' @keywords internal
#' 
#' @noRd
validate_infercnv_allele_obj <- function(infercnv_allele_obj){
  flog.info("validating infercnv_allele_obj")
  
  if(isTRUE(all.equal(colnames(infercnv_allele_obj@allele.data),
                      colnames(infercnv_allele_obj@coverage.data)))){
    if(isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data),
                        rownames(infercnv_allele_obj@coverage.data)))){
      if(isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data),
                          names(infercnv_allele_obj@SNP_info)))){
        if(isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data),
                            rownames(infercnv_allele_obj@gene_order)))){
          if(is.null(infercnv_allele_obj@expr.data) & is.null(infercnv_allele_obj@count.data)){
            return()
          } else{
            if(isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data),
                                rownames(infercnv_allele_obj@expr.data)))){
              if(isTRUE(all.equal(rownames(infercnv_allele_obj@allele.data),
                                  rownames(infercnv_allele_obj@count.data)))){
                return()
              } else{
                flog.error("hmm.... rownames(infercnv_allele_obj@allele.data != 
                           rownames(infercnv_allele_obj@count.data))")
              }
            } else{
              flog.error("hmm.... rownames(infercnv_allele_obj@allele.data != 
                         rownames(infercnv_allele_obj@expr.data))")
            }
          }
          #return()
        } else{
          flog.error("hmm.... rownames(infercnv_allele_obj@allele.data != 
                     rownames(infercnv_allele_obj@gene_order))")
        }
        
      } else{
        flog.error("hmm.... rownames(infercnv_allele_obj@allele.data != 
                   rownames(infercnv_allele_obj@SNP_info))")
      }
      
    } else{
      flog.error("hmm.... rownames(infercnv_allele_obj@allele.data != 
                 rownames(infercnv_allele_obj@coverage.data))")
    }
  } else{
    flog.error("hmm.... colnames(infercnv_allele_obj@allele.data != 
               colnames(infercnv_allele_obj@coverage.data))")
  }
  broken.infercnv_allele_obj = infercnv_allele_obj
  save(broken.infercnv_allele_obj, file="broken.infercnv_allele_obj")
  stop("Problem detected w/ infercnv_allele_obj")
}

#' @title remove_snps()
#'
#' @description infercnv_allele obj accessor method to remove snps from the matrices
#'
#' @param infercnv_allele_obj infercnv_allele object
#' 
#' @param snps_indices_to_remove matrix indices for snps to remove
#'
#' @return infercnv_allele_obj
#'
#' @keywords internal
#' 
#' @noRd
#'
remove_snps <- function(infercnv_allele_obj, snps_indices_to_remove) {
  
  infercnv_allele_obj@expr.data <- infercnv_allele_obj@expr.data[ -1 * snps_indices_to_remove, , drop=FALSE]
  
  infercnv_allele_obj@count.data <- infercnv_allele_obj@count.data[ -1 * snps_indices_to_remove, , drop=FALSE]
  
  infercnv_allele_obj@allele.data <- infercnv_allele_obj@allele.data[ -1 * snps_indices_to_remove, , drop=FALSE]
  
  infercnv_allele_obj@coverage.data <- infercnv_allele_obj@coverage.data[ -1 * snps_indices_to_remove, , drop=FALSE]
  
  infercnv_allele_obj@SNP_info <- infercnv_allele_obj@SNP_info[ -1 * snps_indices_to_remove]
  
  infercnv_allele_obj@gene_order <- infercnv_allele_obj@gene_order[ -1 * snps_indices_to_remove, , drop=FALSE]
  
  validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}

#' @title setAlleleMatrix
#' 
#' @description This function initializes the lesser allele fraction/count determined by at least the minimum
#' number of snps having allele fraction less/more than 0.5 with the minimum of cells supported. 
#' lesser allele fraction will be stored in the slot @expr.data. lesser allele count 
#' will be stored in the slot @count.data
#' 
#' @param infercnv_allele_obj infercnv_allele_object
#' 
#' @param snp_min_coverage the minimum number of threshold that each instance should have in allele/coverage data. 
#' Each instance that is less than this number will be replaced with 0. default = 0 (No action)
#' 
#' @param snp_filter a boolean value whether to do pre-filtering for allele matrix. default = TRUE
#' 
#' @param snp_min.cell a threshold used to filter out non-heterozygous snps. The minimum number of (normal) cells 
#' having both alleles. default = 3
#' 
#' @export
setAlleleMatrix <- function(infercnv_allele_obj,
                            snp_min_coverage = 0,
                            snp_filter = TRUE,
                            snp_min.cell = 3){
  
  if(snp_min_coverage > 1){
    flog.info(sprintf("Replacing allele/coverge instances that are less than %s with zero ...",
                      snp_min_coverage))
    
    infercnv_allele_obj@allele.data[infercnv_allele_obj@coverage.data < snp_min_coverage] <- 0 
    infercnv_allele_obj@coverage.data[infercnv_allele_obj@coverage.data < snp_min_coverage] <- 0
    
    if(mean(rowSums(infercnv_allele_obj@coverage.data) == 0) > 0){
      infercnv_allele_obj <- infercnv:::remove_snps(infercnv_allele_obj, 
                                                    which(rowSums(infercnv_allele_obj@coverage.data) == 0))
    }
  }
  
  flog.info("Estimating allele frequency profiles ...")
  
  allele_matrix <- infercnv_allele_obj@allele.data/infercnv_allele_obj@coverage.data
  allele_matrix[is.na(allele_matrix)] <- 0 # omit no coverage
  allele_matrix[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
  
  if (snp_filter){
    
    if(!is.null(unlist(infercnv_allele_obj@reference_grouped_cell_indices))){
      
      cells_w_ref_allele = rowSums(allele_matrix[,unlist(infercnv_allele_obj@reference_grouped_cell_indices)] != 0 & 
                                     allele_matrix[,unlist(infercnv_allele_obj@reference_grouped_cell_indices)] < 0.5)
      cells_w_alt_allele = rowSums(allele_matrix[,unlist(infercnv_allele_obj@reference_grouped_cell_indices)] != 0 & 
                                     allele_matrix[,unlist(infercnv_allele_obj@reference_grouped_cell_indices)] > 0.5)
      
    } else{
      
      cells_w_ref_allele = rowSums(allele_matrix !=0 & 
                                     allele_matrix < 0.5)
      cells_w_alt_allele = rowSums(allele_matrix !=0 & 
                                     allele_matrix > 0.5)
      
    }
    
    snp_index <- (cells_w_ref_allele >= snp_min.cell & cells_w_alt_allele >= snp_min.cell)
    num_snps_all = sum(snp_index)
    flog.info(sprintf("Number of heterozygous snps used for modeling: %s ...", num_snps_all))
    
    infercnv_allele_obj <- infercnv:::remove_snps(infercnv_allele_obj,which(!snp_index))
    allele_matrix <- allele_matrix[snp_index,]
  }
  
  flog.info("Setting composite lesser allele fraction ...")
  
  lesse_allele_index <- c()
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
    lesse_allele_index <<- c(lesse_allele_index,frac_high)
    if ( frac_high > 0.5) {
      x[x==1] = 0.999
      x[nonzero_val_idx ] = 1 - x[nonzero_val_idx]
    }
    x
  })
  allele_matrix = t(mAF_allele_matrix)
  
  thr_index <- lesse_allele_index > 0.5
  
  lesser.allele.data <- infercnv_allele_obj@allele.data
  lesser.allele.data[thr_index,] <- infercnv_allele_obj@coverage.data[thr_index,] - infercnv_allele_obj@allele.data[thr_index,]
  
  lesser.allele.fraction <- allele_matrix
  infercnv_allele_obj@expr.data <- lesser.allele.fraction
  infercnv_allele_obj@count.data <- lesser.allele.data
  
  # if (smooth_method == 'runmeans') {
  #   
  #   infercnv_allele_obj <- infercnv:::smooth_by_chromosome_runmeans(infercnv_allele_obj,
  #                                                        window_length)
  # } else if (smooth_method == 'pyramidinal') {
  #   
  #   infercnv_allele_obj <- infercnv:::smooth_by_chromosome(infercnv_allele_obj,
  #                                               window_length=window_length,
  #                                               smooth_ends=TRUE)
  # } else if (smooth_method == 'coordinates') {
  #   infercnv_allele_obj <- infercnv:::smooth_by_chromosome_coordinates(infercnv_allele_obj,
  #                                                           window_length=window_length)
  # } else {
  #   stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
  # }
  
  infercnv:::validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}


#' @title setAlleleMatrix_HB -- deprecated
#' 
#' @description This function mimics the way the HB does aiming to initialize the lesser allele fraction/count.
#' lesser allele fraction will be stored in the slot @expr.data. lesser allele count 
#' will be stored in the slot @count.data
#' 
#' @param infercnv_allele_obj infercnv_allele_object
#' 
#' @param snp_min_coverage the minimum number of threshold that each instance should have in allele/coverage data. 
#' Each instance that is less than this number will be replaced with 0. default = 0 (No action)
#' 
#' @param snp_filter a boolean value whether to do pre-filtering for allele matrix. default = TRUE
#' 
#' @param snp_het.deviance.threshold a threshold used to filter 
#' out non-hetergozous snps. default = 0.1
#' 
#' @param snp_min.cell a threshold used to filter out snps that express less than 
#' the minimum number of cells. default = 3
#' 
#' @export
setAlleleMatrix_HB <- function(infercnv_allele_obj,
                               snp_min_coverage = 0,
                               snp_filter = TRUE,
                               snp_het.deviance.threshold = 0.1, 
                               snp_min.cell = 3){
  
  if(snp_min_coverage > 1){
    flog.info(sprintf("Replacing allele/coverge instances that are less than %s with zero ...",
                      snp_min_coverage))
    
    infercnv_allele_obj@allele.data[infercnv_allele_obj@coverage.data < snp_min_coverage] <- 0 
    infercnv_allele_obj@coverage.data[infercnv_allele_obj@coverage.data < snp_min_coverage] <- 0
    
    if(mean(rowSums(infercnv_allele_obj@coverage.data) == 0) > 0){
      infercnv_allele_obj <- infercnv:::remove_snps(infercnv_allele_obj, 
                                                    which(rowSums(infercnv_allele_obj@coverage.data) == 0))
    }
  }
  
  flog.info("Creating in-silico bulk ...")
  allele_bulk <- rowSums(infercnv_allele_obj@allele.data > 0)
  coverage_bulk <- rowSums(infercnv_allele_obj@coverage.data > 0)
  
  if (snp_filter){
    E <- allele_bulk/coverage_bulk
    filter_index_het <- E > snp_het.deviance.threshold & E < 1-snp_het.deviance.threshold
    if(sum(filter_index_het) < 0.01*length(allele_bulk)) {
      flog.info("WARNING! CLONAL DELETION OR LOH POSSIBLE!")
    }
    filter_index_min <- rowSums(infercnv_allele_obj@coverage.data > 0) >= snp_min.cell
    flog.info(sprintf("%s heterozygous SNPs identified ...", 
                      sum(filter_index_het & filter_index_min)))
    
    if(sum(filter_index_het & filter_index_min) < length(allele_bulk)){
      infercnv_allele_obj <- remove_snps(infercnv_allele_obj, 
                                         which(!(filter_index_het & filter_index_min)))
      allele_bulk <- allele_bulk[filter_index_het & filter_index_min]
      coverage_bulk <- coverage_bulk[filter_index_het & filter_index_min]
    }
  }
  
  flog.info("Setting composite lesser allele fraction ...")
  E <- allele_bulk/coverage_bulk
  thr_index <- E > 0.5
  
  lesser.allele.data <- infercnv_allele_obj@allele.data
  lesser.allele.data[thr_index,] <- infercnv_allele_obj@coverage.data[thr_index,] - infercnv_allele_obj@allele.data[thr_index,]
  lesser.allele.fraction <- lesser.allele.data/infercnv_allele_obj@coverage.data
  
  lesser.allele.fraction[is.na(lesser.allele.fraction)] <- 0 # omit no coverage
  lesser.allele.fraction[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
  #lesser.allele.fraction[infercnv_allele_obj@coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
  
  infercnv_allele_obj@expr.data <- lesser.allele.fraction
  infercnv_allele_obj@count.data <- lesser.allele.data
  
  # if (smooth_method == 'runmeans') {
  # 
  #   infercnv_allele_obj <- smooth_by_chromosome_runmeans(infercnv_allele_obj,
  #                                                        window_length)
  # } else if (smooth_method == 'pyramidinal') {
  # 
  #   infercnv_allele_obj <- smooth_by_chromosome(infercnv_allele_obj,
  #                                               window_length=window_length,
  #                                               smooth_ends=TRUE)
  # } else if (smooth_method == 'coordinates') {
  #   infercnv_allele_obj <- smooth_by_chromosome_coordinates(infercnv_allele_obj,
  #                                                           window_length=window_length)
  # } else {
  #   stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
  # }

  validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}


#' @title map2gene
#' 
#' @description This function aims to map snp data back to gene data given the annotation file provided
#' 
#' @param infercnv_allele_obj infercnv_allele based obj
#' 
#' @param gene_order gene annotation file
#' 
#' @export
map2gene <- function(infercnv_allele_obj,
                     gene_order){
  
  flog.info("Mapping snps back to genes ...")
  
  gene_order_gr <- GRanges(gene_order[[C_CHR]],
                           IRanges(as.numeric(as.character(gene_order[[C_START]])),
                                   as.numeric(as.character(gene_order[[C_STOP]]))))
  gene_order_gr$gene <- gene_order %>% rownames()
  
  snp2gene_index <- nearest(infercnv_allele_obj@SNP_info, gene_order_gr)
  infercnv_allele_obj@SNP_info$gene <- gene_order_gr$gene[snp2gene_index]
  
  validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}

## aggregrate into a higher level
# agg_method <- function(data, index, method = c("mean","median"),
#                        cores = 5){
# 
#   method <- match.arg(method)
#   
#   tmp <- mclapply(seq_len(ncol(data)), 
#                   function(x){
#                     tapply(data[,x], 
#                            index, 
#                            method)
#                   },mc.cores = cores)
#   tmp <- do.call(cbind, tmp)
#   colnames(tmp) <- colnames(data)
#   #tmp <- round(tmp)
#   
#   return(tmp)
#   
# }

# collapse_snp2gene
# This function aims to collapse infercnv_allele obj from snp level into gene level
# infercnv_allele_obj infercnv_allele based obj
# gene_annot a gene order annotation
# collapse_method The method used to collapse snp into gene. default = median
# cores A number of cores being used during parallel computing. default = 5
# collapse_snp2gene <- function(infercnv_allele_obj,
#                               gene_annot,
#                               collapse_method = c("mean","median"),
#                               cores = 5){
# 
#   collapse_method <- match.arg(collapse_method)
#   
#   flog.info(sprintf("Using %s method to aggregrate snp into gene level", 
#                     collapse_method))
#   
#   expr.data <- agg_method(data = infercnv_allele_obj@expr.data,
#                           index = infercnv_allele_obj@SNP_info$gene,
#                           method = collapse_method,
#                           cores = cores)
#   
#   count.data <- agg_method(data = infercnv_allele_obj@count.data,
#                            index = infercnv_allele_obj@SNP_info$gene,
#                            method = collapse_method,
#                            cores = cores) %>% round()
#   
#   allele.data <- agg_method(data = infercnv_allele_obj@allele.data,
#                             index = infercnv_allele_obj@SNP_info$gene,
#                             method = collapse_method,
#                             cores = cores) %>% round()
#   
#   coverage.data <- agg_method(data = infercnv_allele_obj@coverage.data,
#                               index = infercnv_allele_obj@SNP_info$gene,
#                               method = collapse_method,
#                               cores = cores) %>% round()
#   
#   gene_order <- gene_annot[rownames(gene_annot) %in% rownames(expr.data),]
#   gene_order_gr <- GenomicRanges::GRanges(gene_order[[C_CHR]],
#                                           IRanges::IRanges(gene_order[[C_START]],
#                                                            gene_order[[C_STOP]]))
#   names(gene_order_gr) <- rownames(gene_order)
#   gene_order_gr <- gene_order_gr %>% sortSeqlevels() %>% sort()
#   
#   expr.data <- expr.data[names(gene_order_gr), ]
#   count.data <- count.data[names(gene_order_gr), ]
#   allele.data <- allele.data[names(gene_order_gr), ]
#   coverage.data <- coverage.data[names(gene_order_gr), ]
#   gene_order <- gene_order[names(gene_order_gr),]
#   gene_order_gr$gene <- names(gene_order_gr)
#   
#   names(gene_order_gr) <- 
#     rownames(expr.data) <- 
#     rownames(count.data)  <-
#     rownames(allele.data) <-
#     rownames(coverage.data) <-
#     rownames(gene_order) <-
#     paste0(gene_order[[C_CHR]], ":",
#            gene_order[[C_START]], ":",
#            gene_order[[C_STOP]])
#   
#   infercnv_allele_obj@expr.data <- expr.data
#   infercnv_allele_obj@count.data <- count.data
#   infercnv_allele_obj@allele.data <- allele.data
#   infercnv_allele_obj@coverage.data <- coverage.data
#   infercnv_allele_obj@SNP_info <- gene_order_gr
#   infercnv_allele_obj@gene_order <- gene_order
#   
#   validate_infercnv_allele_obj(infercnv_allele_obj)
#   
#   return(infercnv_allele_obj)
# }

#' @title collapse_snp2gene
#' 
#' @description This function aims to collapse infercnv_allele obj from snp level into gene level
#' by only leveraging different characteristics
#' 
#' @param infercnv_allele_obj infercnv_allele based obj
#' 
#' @param gene_annot a gene order annotation
#' 
#' @param collapse_method The method used to collapse snp into gene. Take the highest coverage/median/mean
#' value within each gene as its representative. default = "highest"
#' 
#' @export
collapse_snp2gene <- function(infercnv_allele_obj,
                              gene_annot,
                              collapse_method = c("highest","median","mean")){
  
  collapse_method <- match.arg(collapse_method)
  
  flog.info(sprintf("Using %s method to aggregrate snp into gene level", 
                    collapse_method))
  
  melt_data <- melt(infercnv_allele_obj@coverage.data)
  colnames(melt_data) <- c("chrpos","cell","COV")
  
  melt_data$ALT <- melt(infercnv_allele_obj@count.data)[,3]
  melt_data$allele <- melt(infercnv_allele_obj@allele.data)[,3]
  melt_data$AF <- melt(infercnv_allele_obj@expr.data)[,3]
  
  melt_data <- melt_data %>%
    dplyr::filter(COV > 0)
  melt_data$gene <- infercnv_allele_obj@SNP_info[melt_data$chrpos]$gene
  
  if(collapse_method == "median"){
    
    subset_melt_data <- melt_data %>% 
      group_by(cell, gene) %>% 
      mutate(AF = median(AF),
             ALT = median(ALT),
             allele = median(allele),
             COV = median(COV)) %>% 
      ungroup()
    
  } 
  
  if(collapse_method == "mean"){
    
    subset_melt_data <- melt_data %>% 
      group_by(cell, gene) %>% 
      mutate(AF = mean(AF),
             ALT = mean(ALT),
             allele = mean(allele),
             COV = mean(COV)) %>% 
      ungroup()
    
  }
  
  if(collapse_method == "highest"){
    
    subset_melt_data <- melt_data %>% 
      group_by(cell, gene) %>% 
      arrange(desc(COV)) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      ungroup()
    
  }
  
  expr.data <- subset_melt_data %>% select(cell, gene, AF) %>% unique() %>% 
    pivot_wider(names_from = cell,
                values_from = AF,
                values_fill = 0) %>% data.frame()
  rownames(expr.data) <- expr.data$gene
  expr.data[,"gene"] <- NULL
  expr.data <- as.matrix(expr.data)
  
  count.data <- subset_melt_data %>% select(cell, gene, ALT) %>% unique() %>% 
    pivot_wider(names_from = cell,
                values_from = ALT,
                values_fill = 0) %>% data.frame()
  rownames(count.data) <- count.data$gene
  count.data[,"gene"] <- NULL
  count.data <- as.matrix(count.data) %>% round()
  
  allele.data <- subset_melt_data %>% select(cell, gene, allele) %>% unique() %>% 
    pivot_wider(names_from = cell,
                values_from = allele,
                values_fill = 0) %>% data.frame()
  rownames(allele.data) <- allele.data$gene
  allele.data[,"gene"] <- NULL
  allele.data <- as.matrix(allele.data) %>% round()
  
  coverage.data <- subset_melt_data %>% select(cell, gene, COV) %>% unique() %>% 
    pivot_wider(names_from = cell,
                values_from = COV,
                values_fill = 0) %>% data.frame()
  rownames(coverage.data) <- coverage.data$gene
  coverage.data[,"gene"] <- NULL
  coverage.data <- as.matrix(coverage.data) %>% round()
  
  expr.data <- expr.data[unique(infercnv_allele_obj@SNP_info$gene),
                         colnames(infercnv_allele_obj@expr.data)]
  count.data <- count.data[unique(infercnv_allele_obj@SNP_info$gene),
                           colnames(infercnv_allele_obj@expr.data)]
  allele.data <- allele.data[unique(infercnv_allele_obj@SNP_info$gene),
                             colnames(infercnv_allele_obj@expr.data)]
  coverage.data <- coverage.data[unique(infercnv_allele_obj@SNP_info$gene),
                                 colnames(infercnv_allele_obj@expr.data)]
  
  gene_order <- gene_annot[rownames(gene_annot) %in% rownames(expr.data),]
  gene_order_gr <- GenomicRanges::GRanges(gene_order[[C_CHR]],
                                          IRanges::IRanges(gene_order[[C_START]],
                                                           gene_order[[C_STOP]]))
  names(gene_order_gr) <- rownames(gene_order)
  gene_order_gr <- gene_order_gr %>% sortSeqlevels() %>% sort()
  
  expr.data <- expr.data[names(gene_order_gr), ]
  count.data <- count.data[names(gene_order_gr), ]
  allele.data <- allele.data[names(gene_order_gr), ]
  coverage.data <- coverage.data[names(gene_order_gr), ]
  gene_order <- gene_order[names(gene_order_gr),]
  gene_order_gr$gene <- names(gene_order_gr)
  
  names(gene_order_gr) <- 
    rownames(expr.data) <- 
    rownames(count.data)  <-
    rownames(allele.data) <-
    rownames(coverage.data) <-
    rownames(gene_order) <-
    paste0(gene_order[[C_CHR]], ":",
           gene_order[[C_START]], ":",
           gene_order[[C_STOP]])
  
  infercnv_allele_obj@expr.data <- expr.data
  infercnv_allele_obj@count.data <- count.data
  infercnv_allele_obj@allele.data <- allele.data
  infercnv_allele_obj@coverage.data <- coverage.data
  infercnv_allele_obj@SNP_info <- gene_order_gr
  infercnv_allele_obj@gene_order <- gene_order
  
  validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}


#' @title plot_allele
#' 
#' @description plot a summary figure containing allele frequency
#' 
#' @keywords internal
#' 
#' @noRd
plot_allele <- function(infercnv_allele_obj,
                        initialzied_method = c("default","HB"),
                        #expression_mode = F,
                        allele_frequency_mode = F,
                        #HMM = NULL,
                        #use_common_gene = F, 
                        name_to_plot,
                        trend_smK = 31,
                        CELL_POINT_ALPHA = 0.6, dotsize=0.3, colorscheme = "BlueRed"){
  
  initialzied_method <- match.arg(initialzied_method)
  
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
  if(is.null(infercnv_allele_obj@expr.data) | is.null(infercnv_allele_obj@count.data)){
    #flog.error("Please set allele matrix before plotting!")
    flog.info(sprintf("Setting allele matrix using %s way", initialzied_method))
    
    if(initialzied_method == "default"){
      infercnv_allele_obj <- setAlleleMatrix(infercnv_allele_obj)
    } else {
      infercnv_allele_obj <- setAlleleMatrix_HB(infercnv_allele_obj)
    }
    
  }
  
  ############# allele data
  flog.info("Building snp plot ...")
  
  #allele_matrix <- infercnv_allele_obj@allele.data/infercnv_allele_obj@coverage.data # fraction
  #allele_matrix <- infercnv_allele_obj@allele.lesser.data/infercnv_allele_obj@coverage.data
  
  # allele_matrix[is.na(allele_matrix)] <- 0 # omit no coverage
  # allele_matrix[infercnv_allele_obj@allele.data == 0 & infercnv_allele_obj@coverage.data != 0] <- 0.001 # pseudo count for total coverage not 0
  # allele_matrix[infercnv_allele_obj@coverage.data <= 2] <- 0 # filter with the cutoff of 3 coverage
  # 
  # min.cells = 3
  # cells_w_ref_allele = rowSums(allele_matrix != 0 & allele_matrix < 0.5)
  # cells_w_alt_allele = rowSums(allele_matrix != 0 & allele_matrix > 0.5)
  # 
  # allele_matrix = allele_matrix[(cells_w_ref_allele >= min.cells & cells_w_alt_allele >= min.cells), ]
  # num_snps_all = nrow(allele_matrix)
  # flog.info(sprintf("Number of heterozygous snps used for plotting: %s ...", num_snps_all))
  # 
  # flog.info("Setting alt allele fraction to the tumor cell-population minor allele ...")
  # 
  # mAF_allele_matrix = apply(allele_matrix, 1, function(x) {
  #   nonzero_val_idx = which(x>0)
  #   nonzero_vals = x[nonzero_val_idx]
  #   
  #   frac_high = sum(nonzero_vals>0.5)/length(nonzero_vals)
  #   
  #   ## focus allele selection based on the tumor cells only.
  #   tumor_vals = x[unlist(infercnv_allele_obj@observation_grouped_cell_indices)]
  #   tumor_nonzero_vals = tumor_vals[tumor_vals>0]
  #   if (length(tumor_nonzero_vals) > 0) {
  #     frac_high = sum(tumor_nonzero_vals>0.5)/length(tumor_nonzero_vals)
  #   }
  #   
  #   if ( frac_high > 0.5) {
  #     x[x==1] = 0.999
  #     x[nonzero_val_idx ] = 1 - x[nonzero_val_idx]
  #   }
  #   x
  # })
  # 
  # allele_matrix = t(mAF_allele_matrix)
  # 
  allele_matrix = infercnv_allele_obj@expr.data
  
  if (!is.null(infercnv_allele_obj@tumor_subclusters)) {
    ## define cell ordering.
    flog.info("Exracting clustering info ...")
    #ordered_cells <- colnames(allele_matrix)[unlist(infercnv_allele_obj@tumor_subclusters$subclusters)]
    tmp_sub <- lapply(infercnv_allele_obj@tumor_subclusters$subclusters, rev)
    ordered_cells <- colnames(allele_matrix)[unlist(tmp_sub)]
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
  #gexdata <- NULL
  #chr_maxpos_gex <- NULL
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
  
  make_plots = function(alleledataToPlot, 
                        #gexdataToPlot = NULL, 
                        chr_maxpos_snp = chr_maxpos_snp, 
                        #chr_maxpos_gex = chr_maxpos_gex,
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
      scale_color_manual(values = color_cell) + ylim(c(0,1))
    
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
    
    if(allele_frequency_mode){

      cell_count <- alleledatamelt %>% select(chr, cell, sample_type) %>% unique() %>%
        group_by(chr, sample_type) %>% summarise(cell_count = n(), .groups = 'drop') %>%
        select(cell_count)
      bar_data <- alleledatamelt %>% select(chrpos, chr, sample_type) %>% unique() %>%
        group_by(chr, sample_type) %>% summarise(snp_count = n(), .groups = 'drop') %>%
        mutate(cell_count, "SNPs coverage per cell" = snp_count/cell_count)

      bar_plot <- bar_data %>% ggplot(aes(x = sample_type, y = `SNPs coverage per cell`, fill = sample_type)) +
        geom_bar(stat='identity') +
        facet_grid (~chr, scales = 'free_x', space = 'fixed') +

        labs(title="Fraction of heterozygous snps used") +
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
        scale_fill_manual(values = color_cell)

    }
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
      
      if(allele_frequency_mode){
        
        ratio_normal_cells = max(0.5, num_normal_cells/(num_normal_cells + num_malignant_cells))
        pg = plot_grid(normal_snps_plot,
                       bar_plot,
                       allele_freq_plot_w_trendlines, 
                       malignant_snps_plot, 
                       ncol=1, align='v', 
                       rel_heights=c(ratio_normal_cells, 
                                     0.25,
                                     0.25,
                                     1-ratio_normal_cells))
        
      } else{
        
        ratio_normal_cells = max(0.25, num_normal_cells/(num_normal_cells + num_malignant_cells))
        pg = plot_grid(normal_snps_plot, 
                       allele_freq_plot_w_trendlines, 
                       malignant_snps_plot, 
                       ncol=1, align='v', 
                       rel_heights=c(ratio_normal_cells, 0.25, 1-ratio_normal_cells)) 
      }
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
                  #gexdataToPlot = gexdata,
                  chr_maxpos_snp = chr_maxpos_snp, 
                  #chr_maxpos_gex = chr_maxpos_gex, 
                  infercnv_allele_obj = infercnv_allele_obj)
  
  ggsave (name_to_plot, p, width = 13.33, height = 7.5, units = 'in', dpi = 300)
  flog.info("Done!")
}
