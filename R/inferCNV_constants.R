#!/usr/bin/env Rscript



C_CHR <- "chr"
C_START <- "start"
C_STOP <- "stop"
C_HCLUST_METHODS <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
C_OUTPUT_FORMAT <- c("pdf", "png")


## also including some globals:
infercnv.env <- new.env()
infercnv.env$GLOBAL_NUM_THREADS <- 1  # default is single-threaded.


#' @importFrom grDevices col2rgb colorRampPalette dev.off pdf png rgb
#' @importFrom graphics abline axis boxplot hist image layout lines mtext par plot points rect text title legend
#' @importFrom stats as.dendrogram as.dist as.hclust cutree density dist filter median order.dendrogram quantile reorder sd complete.cases cor t.test p.adjust predict rnorm runif smooth.spline var wilcox.test dnorm ecdf ks.test lm nls pnorm qgamma qnorm rbinom rchisq rgamma rlnorm rnbinom rpois shapiro.test update
#' @importFrom utils flush.console read.table write.table tail read.csv head
#' @import futile.logger
#' @importFrom methods setClass new is
#' @importFrom gplots bluered
#' @importFrom ape write.tree as.phylo read.tree drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @importFrom fastcluster hclust
#' @importFrom parallelDist parallelDist
#' @import RColorBrewer
#' @importFrom Matrix Matrix rowMeans colSums sparseMatrix
#' @importFrom dplyr %>% count
#' @import fitdistrplus
#' @import foreach
#' @import doParallel
#' @import future
#' @import coda
#' @import ggplot2
#' @import argparse
#' @importFrom edgeR estimateDisp
#' @importFrom caTools runmean
#' @importFrom coin oneway_test pvalue
#' @importFrom digest digest
#' @importFrom RANN nn2
#' @importFrom igraph graph_from_adjacency_matrix cluster_leiden
#' @importFrom reshape2 melt
#' @importFrom rjags jags.model coda.samples
#' @importFrom BiocGenerics counts t mean
#' @importFrom SummarizedExperiment start colData rowData assays assays<- rowData<- colData<- 
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom tidyr gather
#' @importFrom parallel detectCores
#' @importFrom gridExtra ttheme_default tableGrob gtable_combine marrangeGrob 
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges `%over%`
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom S4Vectors queryHits subjectHits 


NULL

