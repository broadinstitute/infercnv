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
  min_separation = 1,
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

\item{min_separation}{minimum difference of the two k-means centres of log CNV score in pass 1 (default: 1, i.e. an
e-fold difference); below it no distinct normal population is assumed and the function stops}

\item{...}{further arguments to run() for both passes (e.g. num_threads, k_nn, leiden_resolution)}
}
\value{
list(infercnv_obj = object with the new reference groups, cells = data.frame per cell: cell, group,
subcluster, class (reference / normal / malignant / unresolved), cnv_score, cnv_cor)
}
\description{
Finds normal (copy-number neutral) cells to use as the reference when none are known, and classifies the remaining
cells as malignant or normal, in two passes of inferCNV up to the subclustering step (no HMM). The provisional reference
comes from normal_groups if given; otherwise from subclusters expressing markers of non-malignant lineages (immune,
endothelial, fibroblast), as SCEVAN does (De Falco et al. 2023, Nat Commun 14:1074); failing both, from the subclusters
with the flattest copy-number profiles, the CopyKAT criterion (Gao et al. 2021, Nat Biotechnol 39:599), which is
unreliable when most cells are tumour cells (the tumour then sets the baseline). Pass 2 reruns inferCNV with that
reference and classifies each subcluster by CNV score and CNV correlation as in Tirosh et al. 2016 (Science 352:189) and
Puram et al. 2017 (Cell 171:1611), with thresholds from random groups of reference cells of the same size. A warning is
given if most cells called malignant express non-malignant markers (an inverted reference). Marker panels assume the
tumour is not of an immune, endothelial or fibroblast lineage.
}
