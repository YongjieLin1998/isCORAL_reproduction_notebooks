\name{CalculateGeneOmegaSquared}
\alias{CalculateGeneOmegaSquared}
\title{Calculate gene-wise omega-squared effect sizes}
\usage{
CalculateGeneOmegaSquared(
  expression_matrix,
  groups,
  chunk_size = 2000,
  verbose = FALSE
)
}
\arguments{
  \item{expression_matrix}{Gene-by-cell expression matrix. Sparse matrices are supported.}
  \item{groups}{Cell-group assignments aligned with the matrix columns. Missing assignments are ignored.}
  \item{chunk_size}{Number of genes processed in each memory-bounded chunk.}
  \item{verbose}{Whether to report degenerate grouping conditions.}
}
\value{A named numeric vector containing one omega-squared value per gene.}
\description{
Computes a one-way, gene-wise omega-squared effect size for the supplied cell
groups while processing genes in chunks.
}

