% Please edit documentation in R/inferCNV_auto_reference.R
\name{infercnv_auto_reference}
\alias{infercnv_auto_reference}
\title{infercnv_auto_reference}
\usage{
infercnv_auto_reference(
  infercnv_obj,
  out_dir,
  cutoff = 1,
  normal_groups = NULL,
  normal_markers = .default_normal_markers,
  marker_fold = 4,
  score_quantile = 0.99,
  min_ref_cells = 20,
  reference_method = c("auto", "markers", "self_centered"),
  flat_tolerance = 1,
  min_flat_cells = 10,
  flat_resolution = 1,
  flat_n_features = 2000,
  ...
)
}
\arguments{
\item{infercnv_obj}{infercnv object (references, if any, are ignored and treated as observations)}

\item{out_dir}{directory for the two passes (subdirectories auto_reference_pass1 and _pass2)}

\item{cutoff}{expression cutoff passed to run() (1 for Smart-seq2, 0.1 for 10x)}

\item{normal_groups}{annotation groups known to be normal (optional; used as the provisional reference)}

\item{normal_markers}{named list of marker gene panels of non-malignant lineages (default: immune, endothelial,
fibroblast), or NULL to skip the marker step}

\item{marker_fold}{minimum fold difference of mean panel expression between marker-positive and other subclusters
(default: 4)}

\item{score_quantile}{quantile of the reference-group null used as the threshold (default: 0.99)}

\item{min_ref_cells}{minimum number of provisional reference cells (default: 20)}

\item{reference_method}{how the provisional reference is found when normal_groups is not given: "auto" (markers
if enough marker-positive cells, otherwise self_centered), "markers" or "self_centered"}

\item{flat_tolerance}{self_centered: subclusters whose flatness score is at most flat_tolerance times that of the
flattest subcluster join the provisional reference (default: 1, the flattest only)}

\item{min_flat_cells}{self_centered: minimum cells for a cluster to be a candidate (default: 10)}

\item{flat_resolution}{self_centered: clustering resolution of the expression clusters (default: 1)}

\item{flat_n_features}{self_centered: number of variable genes for the expression clustering (default: 2000)}

\item{...}{further arguments to run() for both passes (e.g. num_threads, k_nn, leiden_resolution)}
}
\value{
list(infercnv_obj = object with the new reference groups, cells = data.frame per cell: cell, group,
subcluster, class (reference / normal / malignant / unresolved), cnv_score, cnv_cor, marker_positive,
pass1_subcluster, expr_cluster), reference_method ("normal_groups", "markers" or "self_centered"), flatness =
data.frame per expression cluster: cluster, cells, flat_score, ratio_to_flattest, marker_positive, selected)
}
\description{
Finds normal (copy-number neutral) cells to use as the reference when none are known, and classifies the remaining
cells as malignant or normal, in two passes of inferCNV up to the subclustering step (no HMM). The provisional reference
comes from normal_groups if given; otherwise from subclusters expressing markers of non-malignant lineages (immune,
endothelial, fibroblast), as SCEVAN does (De Falco et al. 2023, Nat Commun 14:1074); failing both (or with
reference_method = "self_centered"), from the expression cluster (standard Seurat clustering, i.e. cell types) with
the flattest self-centred profile, the principle of
CopyKAT's baseline search (Gao et al. 2021, Nat Biotechnol 39:599): each subcluster's mean log2(CPM + 1) profile is
centred on its own genome-wide mean (not on the mean of all cells, which a tumour majority would dominate), smoothed
along each chromosome and scored by its mean square. Subclusters within flat_tolerance times the flattest score are
added (a dial for exploring which cell types look normal); all scores are written to auto_reference.flatness.tsv. The
score is not CNV-specific, so in tumours with few CNVs a tumour subcluster can be the flattest. Pass 2 reruns inferCNV with that
reference and classifies each subcluster by CNV score and CNV correlation as in Tirosh et al. 2016 (Science 352:189) and
Puram et al. 2017 (Cell 171:1611), with thresholds from random groups of reference cells of the same size. A warning is
given if most cells called malignant express non-malignant markers (an inverted reference). Marker panels assume the
tumour is not of an immune, endothelial or fibroblast lineage.
}
