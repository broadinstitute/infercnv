Visit project [wiki](https://github.com/broadinstitute/inferCNV/wiki) for InferCNV documentation.

## Run allele_based method:

```
sample_id <- "XXX"

## allele data
alt.matrix <- read.table(path to your alt matrix)
tot.matrix <- read.table(path to your tot matrix)
                                
## expression data
express_data <- read.table(path to your expression data)
annot <- read.table(path to your cell annotation)

## processing
common_cell <- intersect(colnames(express_data),
                         colnames(alt.matrix))

library(infercnv)

save_path <- path_to_your_dir

## create a infercnv obj
infercnv_obj_example = infercnv::CreateInfercnvObject(raw_counts_matrix = express_data[, common_cell],
                                                      raw_allele_matrix = alt.matrix[, common_cell],
                                                      raw_coverage_matrix = tot.matrix[, common_cell],
                                                      gene_order_file = path to your gene annotation,
                                                      annotations_file = annot[rownames(annot) %in% common_cell,,drop = F],
                                                      ref_group_names=c("Microglia/Macrophage"),
                                                      max_cells_per_group = 100,
                                                      snp_split_by = "::")
                                                      
## set allele matrix                                                      
infercnv_obj_example_test <- infercnv::setAlleleMatrix(infercnv_obj_example@.allele)

## do clustering for allele matrix. Here we use leiden method for example.
infercnv_obj_example_cluster_leiden <- infercnv:::define_signif_tumor_subclusters(infercnv_obj_example_test,
                                                                                  cluster_by_groups = F,
                                                                                  partition_method='leiden')

## HMM model using whole sample mode
infercnv_obj_example_HMM_whole_samples <- infercnv:::allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj_example_cluster_leiden) 

## Prediction plot for HMM
infercnv::plot_cnv(infercnv_obj_example_HMM_whole_samples,
                   out_dir=paste0(save_path,sample_id),
                   obs_title="Observations (Cells)",
                   ref_title="References (Cells)",
                   k_obs_groups = 1,
                   cluster_by_groups=F,
                   cluster_references=F,
                   x.center = -0.5,
                   x.range=c(-1,0),
                   color_safe_pal=F,
                   title="Deletion/LOH boundaries (sample mode)",
                   output_filename="sample_mode",
                   output_format="png",
                   png_res=300,
                   dynamic_resize=0,
                   custom_color_pal = color.palette(c("darkblue", "white")))

## subcluster mode
infercnv_obj_example_HMM_subcluster <- infercnv:::allele_HMM_predict_CNV_via_HMM_on_tumor_subclusters(infercnv_obj_example_cluster_leiden)

## prediction plot for HMM
infercnv::plot_cnv(infercnv_obj_example_HMM_subcluster,
                   out_dir=paste0(save_path,sample_id),
                   obs_title="Observations (Cells)",
                   ref_title="References (Cells)",
                   k_obs_groups = 1,
                   cluster_by_groups=F,
                   cluster_references=F,
                   x.center = -0.5,
                   x.range=c(-1,0),
                   color_safe_pal=F,
                   title="Deletion/LOH boundaries (subcluster mode)",
                   output_filename="subcluster_mode",
                   output_format="png",
                   png_res=300,
                   dynamic_resize=0,
                   custom_color_pal = color.palette(c("darkblue", "white")))

## make a snp plot
infercnv:::plot_allele(infercnv_obj_example_cluster_leiden,
                       name_to_plot = plot name with file format)

```

