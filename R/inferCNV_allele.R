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


#' @title setAlleleMatrix_HB
#' 
#' @description This function mimics the way the HB does aiming to initialize the lesser allele fraction/count.
#' lesser allele fraction will be stored in the slot @expr.data. lesser allele count 
#' will be stored in the slot @count.data
#' 
#' @param infercnv_allele_obj infercnv_allele_object
#' 
#' @param snp_min_coverage the minimum number of counts that the snp express in coverage data. default = 2
#' 
#' @param snp_filter a boolean value whether to do pre-filtering for allele matrix. default = TRUE
#' 
#' @param snp_het.deviance.threshold a threshold used to filter 
#' out non-hetergozous snps. default = 0.1
#' 
#' @param snp_min.cell a threshold used to filter out snps that express less than 
#' the minimum number of cells. default = 3
#' 
#' @param smooth_method a method used for smoothing allele fraction default = "runmeans"
#' 
#' @param window_length smoothing window size. default = 101
#' 
#' @export
setAlleleMatrix_HB <- function(infercnv_allele_obj,
                               snp_min_coverage = 2,
                               snp_filter = TRUE,
                               snp_het.deviance.threshold = 0.1, snp_min.cell = 3,
                               smooth_method = "runmeans", window_length = 101){
  
  flog.info("Removing snps that have low coverage ...")
  
  infercnv_allele_obj@allele.data[infercnv_allele_obj@coverage.data <= snp_min_coverage] <- 0 
  infercnv_allele_obj@coverage.data[infercnv_allele_obj@coverage.data <= snp_min_coverage] <- 0  
  
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
  
  if (smooth_method == 'runmeans') {

    infercnv_allele_obj <- smooth_by_chromosome_runmeans(infercnv_allele_obj,
                                                         window_length)
  } else if (smooth_method == 'pyramidinal') {

    infercnv_allele_obj <- smooth_by_chromosome(infercnv_allele_obj,
                                                window_length=window_length,
                                                smooth_ends=TRUE)
  } else if (smooth_method == 'coordinates') {
    infercnv_allele_obj <- smooth_by_chromosome_coordinates(infercnv_allele_obj,
                                                            window_length=window_length)
  } else {
    stop(sprintf("Error, don't recognize smoothing method: %s", smooth_method))
  }

  validate_infercnv_allele_obj(infercnv_allele_obj)
  
  return(infercnv_allele_obj)
}


#' @title map2gene
#' 
#' @description This function aims to map snp data back to gene data given the annotation file
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
