\name{make_PCs_gds}
\alias{make_PCs_gds}
\title{Run Principal Component Analysis (PCA) using the Genomic Data Structure (GDS).}
\usage{
  make_PCs_gds(gdsFile, n.core, n.top = 0)
}
\arguments{
  \item{gdsFile}{A Genomic Data Structure file describing
  your study.}

  \item{n.top}{Number of top principal components to
  return. Defaults to returning all components (i.e. # of
  samples).}

  \item{n.core}{Distrubute computation across N cores.}
}
\value{
  A matrix of Principal Components of dimension (# of
  samples) x (n.top).
}
\description{
  An efficient method for computing Principal Components
  using the Genomic Data Structure (GDS). This is a
  convenience wrapper for functions from the SNPRelate
  package.
}
\references{
  library(SNPRelate)
}
\keyword{GRM}
\keyword{PCA,}
\keyword{covariate,}

