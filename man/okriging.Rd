\name{okriging}
\alias{okriging}
\title{Run omic kriging on a set of correlation matrices and a given phenotype.}
\usage{
  okriging(idtest, idtrain = NULL, corlist, H2vec, pheno,
    phenoname, Xcova = NULL)
}
\arguments{
  \item{corlist}{A list of correlation matrices used in
  Kriging. rownames and colnames of cor should be IID list
  and include idtest and idtrain.}

  \item{H2vec}{has weights for each RM relatednes matrix}

  \item{idtest}{A vector of sample IDs which constitute the
  test set.}

  \item{idtrain}{A vector of sample IDs which constitute
  the training set.}

  \item{pheno}{A data frame with rownames set as sample IDs
  and a column containing phenotype data.}

  \item{phenoname}{The name of the column in pheno which
  contains phenotype data to test.}

  \item{Xcova}{TODO:: What is this?}
}
\value{
  A dataframe with three columns: sample ID, observed
  phenotype Ytest, and predicted phenotype Ypred
}
\description{
  Universal kriging formula: lambda' = ( c + X m )' iSig m'
  = ( x - X' iSig c )' ( X' iSig X )^-1 m' = ( t(x) - c'
  iSig X ) ( X' iSig X )^-1 lambda' = (c' + m' X) iSig x:
  #covariates x ntest X: ntrain x #cov c: ntrain x ntest
}
\references{
  Cressie 1993 Statistics for Spatial Data p.154
}
\keyword{prediction}

