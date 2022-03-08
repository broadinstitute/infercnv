Visit project [wiki](https://github.com/broadinstitute/inferCNV/wiki) for InferCNV documentation.

## Run allele_based method to infer candidate boundaries (HMM):

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

## create a infercnv obj (only select 50 cells for each group)
infercnv_obj_example = infercnv::CreateInfercnvObject(raw_counts_matrix = express_data[, common_cell],
                                                      raw_allele_matrix = alt.matrix[, common_cell],
                                                      raw_coverage_matrix = tot.matrix[, common_cell],
                                                      gene_order_file = path to your gene annotation,
                                                      annotations_file = annot[rownames(annot) %in% common_cell,,drop = F],
                                                      ref_group_names=c("Microglia/Macrophage"),
                                                      max_cells_per_group = 50,
                                                      snp_split_by = "::")
                                                      
## set allele matrix                                                      
infercnv_obj_example_allele <- infercnv::setAlleleMatrix(infercnv_obj_example@.allele)

## do clustering for allele matrix. Here we use leiden method for example.
infercnv_obj_example_cluster_leiden <- infercnv:::define_signif_tumor_subclusters(infercnv_obj_example_allele,
                                                                                  cluster_by_groups = F,
                                                                                  partition_method='leiden')

## do sample mode HMM for each group                                                                                   
infercnv_obj_example_HMM_whole_samples_list <- infercnv:::allele_HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj_example_cluster_leiden)
                                                                                  
## Prediction plot for sample mode HMM
infercnv::plot_cnv(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj,
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
infercnv::plot_cnv(infercnv_obj_example_HMM_subcluster$infercnv_allele_obj,
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

## Run allele-based only bayesian model: -- under development (still many rep codes out there)
```
mcmc_sample_mode_tumor <- infercnv:::inferCNVAlleleBayesNet(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj,
                                                            infercnv_obj_example_HMM_whole_samples_list$HMM_output,                                                                                                    rep(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj@observation_grouped_cell_indices,
                                                            length(infercnv_obj_example_HMM_whole_samples_list$HMM_output)))
                                                                    
mcmc_sample_mode_normal <- infercnv:::inferCNVAlleleBayesNet(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj,
                                                             infercnv_obj_example_HMM_whole_samples_list$HMM_output,
                                                            rep(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj@reference_grouped_cell_indices,
                                                             length(infercnv_obj_example_HMM_whole_samples_list$HMM_output)))
                                                             
## you can interrgrote the results using mcmc_sample_mode_tumor/mcmc_sample_mode_normal@posterior_prob
```

## Run combined method leveraging both gene and allele -- under development (not formalize yet)
```
## have to run regular infercnv step first but with limited steps to get intermediate results
infercnv_obj_example_gene_obj <- infercnv::run(
  infercnv_obj_example,
  cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=paste0("test_result/com_test_gene_", sample_id),
  plot_chr_scale = F,
  cluster_by_groups=TRUE,
  plot_steps=T,
  denoise=TRUE,
  HMM=T,
  analysis_mode = "samples",
  HMM_type = "i3",
  no_prelim_plot=F,
  up_to_step = 17)

## get gene_based hmm obj ("test_result/com_test_gene_" is the prefix of the directory where you run infercnv)  
hmm_obj <- readRDS(paste0("test_result/com_test_gene_", sample_id,
                          "/17_HMM_pred",
                          "HMM",infercnv_obj_example_gene_obj@options$HMM_type,".hmm_mode-",infercnv_obj_example_gene_obj@options$analysis_mode,
                          ".infercnv_obj"))
                          
## create mcmc gene-based obj (though we assign "BUGS_Mixture_Model_i3" as to the bug file,
but it's just for creating step - we don't use it)
MCMC_infercnv_obj_example <- new("MCMC_inferCNV")
MCMC_infercnv_obj_example <- infercnv:::initializeObject(MCMC_infercnv_obj_example, 
                                                          list("file_dir" = paste0("test_result/com_test_gene_", sample_id),
                                                               "resume_file_token" = paste0(paste0("HMM",infercnv_obj_example_gene_obj@options$HMM_type), ".hmm_mode-", infercnv_obj_example_gene_obj@options$analysis_mode),
                                                               "model_file" = system.file("BUGS_Mixture_Model_i3",package = "infercnv"),
                                                               "HMM_type" = infercnv_obj_example_gene_obj@options$HMM_type,
                                                               "quietly" = F), 
                                                          infercnv_obj_example_gene_obj)

## get mean and sd for each state                                                          
MCMC_infercnv_obj_example <- infercnv:::MeanSD(obj = MCMC_infercnv_obj_example, 
                                                HMM_states = hmm_obj@expr.data, 
                                                infercnv_obj = infercnv_obj_example_gene_obj)
                                                
## run mcmc (note: just for sample mode and only for tumor cells for testing)
## it will return a list containing the posterior prob of each disjoin region (allele/gene) having event across chrs
infercnv_obj_example_samples_combined <- infercnv:::inferCNVCombinedBayesNet(infercnv_obj_example_HMM_whole_samples_list$infercnv_allele_obj,
                                                                             infercnv_obj_example_HMM_whole_samples_list$HMM_output,
                                                                             MCMC_infercnv_obj_example,
                                                                             cores = 8)

```