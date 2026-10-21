\name{FindMetaclone}
\alias{FindMetaclone}
\title{Identify Metaclone communities in the CORAL clone space}
\usage{
FindMetaclone(
  seurat_obj,
  resolution = 0.7,
  reduction = "CORAL_clone_embed",
  dims = 1:30,
  graph_prefix = "pcoral",
  cluster_col_name = "metaclone",
  algorithm = 4,
  random_seed = 42,
  verbose = TRUE
)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing the requested clone-space reduction.}
  \item{resolution}{Community-detection resolution passed to \code{Seurat::FindClusters()}.}
  \item{reduction}{Name of the CORAL clone-space reduction.}
  \item{dims}{Dimensions used to construct the nearest-neighbor graph.}
  \item{graph_prefix}{Prefix for the generated nearest-neighbor and shared-nearest-neighbor graphs.}
  \item{cluster_col_name}{Metadata column receiving the Metaclone labels.}
  \item{algorithm}{Clustering algorithm passed to \code{Seurat::FindClusters()}; \code{4} requests Leiden and \code{1} requests Louvain.}
  \item{random_seed}{Random seed for community detection.}
  \item{verbose}{Whether to report graph construction and clustering progress.}
}
\value{The modified Seurat object with Metaclone labels and neighbor graphs.}
\description{
Builds a shared-nearest-neighbor graph in the CORAL clone embedding and
identifies Metaclone communities with Leiden or Louvain clustering. If Leiden
is requested but no supported Leiden backend is available, the function warns
and falls back to Louvain.
}

