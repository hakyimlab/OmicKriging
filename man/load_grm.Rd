\name{load_grm}
\alias{load_grm}
\title{Load genetic relatedness matrix from file.}
\usage{
  load_grm(grmDataFile)
}
\arguments{
  \item{grmDataFile}{Location of the .Rdata file.}
}
\value{
  A genetic correlation matrix with colnames and rownames
  set to sample IDs.
}
\description{
  Loads a genetic relatedness matrix (GRM) from an .Rdata
  file produced by save() in the \code{\link{make_grm}}
  function. This can be used to quickly load a GRM instead
  of recomputing it from PLINK binary files
}
\keyword{input}

