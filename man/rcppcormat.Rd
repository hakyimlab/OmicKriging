\name{rcppcormat}
\alias{rcppcormat}
\title{Compute a correlation matrix.}
\usage{
  rcppcormat(snpmat)
}
\arguments{
  \item{snpmat}{A centered SNP matrix containing fields of
  type 'double'.}
}
\value{
  A correlation matrix.
}
\description{
  This function computes a correlation matrix using a cross
  product. It is implemented in C++ for performance. This
  function is compiled and called from R using the Rcpp
  package.
}
\keyword{c++,}
\keyword{correlation}
\keyword{matrix}
\keyword{performance,}

