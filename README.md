Visit project [wiki](https://github.com/broadinstitute/inferCNV/wiki) for InferCNV documentation.

Run allele_based method:

`sample_id <- "MGH60" `
sample_id <- "MGH60"
### allele data
alt.matrix <- read.table(paste0("data/allele_data/",sample_id,
                                ".alt.dense.matrix"))
tot.matrix <- read.table(paste0("data/allele_data/",sample_id,
                                ".tot.dense.matrix"))
`
###

### expression data
express_data <- read.table(connection <- gzfile("data/gex_data/oligodendro.ALL.rsem.gene_counts.matrix.gz", 'rt'), 
                           header=TRUE, row.names=1, check.names=FALSE)
annot <- read.table("data/gex_data/oligodendro.ALL.cell_annots.txt", 
                    header=FALSE, row.names=1, sep="\t", 
                    stringsAsFactors=FALSE, colClasses = c('character', 'character'))
###

### processing
common_cell <- intersect(colnames(express_data),
                         colnames(alt.matrix))
###
