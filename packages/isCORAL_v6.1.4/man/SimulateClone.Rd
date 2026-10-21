\name{SimulateClone}
\alias{SimulateClone}
\alias{SimulateMetaClone}
\title{Generate hierarchical clone assignments}
\usage{
SimulateClone(seurat_obj, bin_size, k_values, hclust_method = "ward.D2")

SimulateMetaClone(...)
}
\arguments{
  \item{seurat_obj}{A Seurat object containing isCORAL scan results.}
  \item{bin_size}{The bin size identifying an entry in \code{CORAL_scan_results}.}
  \item{k_values}{Numeric vector of cluster counts used to cut the hierarchical tree.}
  \item{hclust_method}{Agglomeration method passed to \code{stats::hclust()}.}
  \item{...}{Arguments forwarded to \code{SimulateClone()}.}
}
\value{
The modified Seurat object. Each requested \code{k} adds a metadata column
named \code{metaclone_k<k>_bs<bin_size>}.
}
\description{
Converts the isCORAL consensus similarity matrix into hierarchical clone
assignments. This is the retained hierarchical-clustering method.
}
\details{
\code{SimulateMetaClone()} is a deprecated optional alias that calls
\code{SimulateClone()}. Embedding-based Metaclone communities are inferred
with \code{FindMetaclone()}.
}

