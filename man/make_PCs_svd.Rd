\name{make_PCs_svd}
\alias{make_PCs_svd}
\title{Run Principal Component Analysis (PCA) using base R svd() function.}
\usage{
  make_PCs_svd(X, n.top = 2)
}
\arguments{
  \item{X}{A correlation matrix.}

  \item{n.top}{Number of top principal compenents to
  return}
}
\value{
  A matrix of Principal Components of dimension (# of
  samples) x (n.top). As expected, eigenvectors are ordered
  by eigenvalue. Rownames are given as sample IDs.
}
\description{
  A simple wrapper around the base R svd() function which
  returns the top N eigenvectors of a matrix. Use this
  function to generate covariates for use with the
  \code{\link{okriging}} or
  \code{\link{krigr_cross_validation}} functions.
}
\keyword{GRM}
\keyword{PCA,}
\keyword{covariate,}

