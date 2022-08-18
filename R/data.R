#' Generated SmartSeq2 expression data with 10 normal cells and 10 tumor cells.
#' This is only to demonstrate how to use methods, not actual data to be used in an analysis.
#'
#' @format A data frame with 1150 rows (genes) and 20 columns (cells)
#'
#'
"infercnv_data_example"

#' Generated SmartSeq2 snp data (alternative) with 10 normal cells and 10 tumor cells.
#' This is only to demonstrate how to use methods, not actual data to be used in an analysis.
#'
#' @format A data frame with 97 rows (snps) and 20 columns (cells)
#'
#'
"infercnv_allele_alt_example"

#' Generated SmartSeq2 snp data (coverage) with 10 normal cells and 10 tumor cells.
#' This is only to demonstrate how to use methods, not actual data to be used in an analysis.
#'
#' @format A data frame with 97 rows (snps) and 20 columns (cells)
#'
#'
"infercnv_allele_tot_example"

#' Generated classification for 10 normal cells and 10 tumor cells.
#'
#' @format A data frame with 20 rows (cells) and 1 columns (classification)
#'
#'
"infercnv_annots_example"

#' Downsampled gene coordinates file from GrCh37
#'
#' @format A data frame with 10338 rows (genes) and 3 columns (chr, start, end)
#'
#'
"infercnv_genes_example"

#' infercnv object result of the processing of run() in the example, to be used for other examples.
#'
#' @format An infercnv object
#'
#'
"infercnv_object_example"

#' infercnv_allele object result of the processing of run() in the example, to be used for other examples.
#'
#' @format An infercnv_allele object
#'
#'
"infercnv_object_allele_example"

#' infercnv_allele object result of the processing of run() in the example, to be used for other examples.
#'
#' @format An infercnv_allele object (gene level)
#'
#'
"infercnv_object_allele_gene_example"

#' infercnv object result of the processing of run() in the HMM example, to be used for other examples.
#'
#' @format An infercnv object containing HMM predictions
#'
#'
"HMM_states"

#' infercnv_allele object result of the processing of run() in the HMM example, to be used for other examples.
#'
#' @format An infercnv_allele object containing HMM predictions
#'
#'
"HMM_allele_states"

#' infercnv_allele object result of the processing of run() in the HMM example, to be used for other examples.
#'
#' @format An infercnv_allele object containing HMM predictions (gene level)
#'
#'
"HMM_allele_gene_states"

#' infercnv object result of the processing of inferCNVBayesNet in the example, to be used for other examples.
#'
#' @format An infercnv object containing posterior probability of CNV states
#'
#'
"mcmc_obj"
