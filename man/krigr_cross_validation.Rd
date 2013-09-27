\name{krigr_cross_validation}
\alias{krigr_cross_validation}
\title{Run n-fold cross validated kriging across n-cores.}
\usage{
  krigr_cross_validation(n.cores, corlist, pheno.df,
    pheno.name)
}
\arguments{
  \item{n.core}{The number of cores used and test sets the
  data will be split into.}

  \item{corlist}{A list of correlation matrices. Colnames
  and rownames must be sample IDs, and must match betweenn
  matrices.}

  \item{pheno.df}{A data frame of sample IDs and their
  corresponding phenotypes. The data frame must have
  rownames() set to sample IDs.}

  \item{pheno.name}{The name of the column in the pheno.df
  which is the phenotype being predicted.}
}
\value{
  A data frame with a column for each: sample IDs, observed
  phenotypes, and predicted phenotypes.
}
\description{
  For now, this is a placeholder function. What is really
  needed is a function which can run across multiple cores
  with an arbitrary number of test/train groups, i.e.
  decouple the n-core and n-fold variables. At its core,
  this function calls \code{\link{okriging}}.
}
\keyword{cross}
\keyword{prediction,}
\keyword{validation}

