\name{krigr_cross_validation}
\alias{krigr_cross_validation}
\title{Multithreaded cross validation routine for Omic Kriging.}
\usage{
  krigr_cross_validation(corlist, pheno.df, pheno.name,
    Xcovamat = NULL, H2vec, nfold = 10, ncore = "all",
    AUC = FALSE, ...)
}
\arguments{
  \item{corlist}{A list of correlation matrices used in
  Kriging. rownames and colnames of cor should be IID list
  and include idtest and idtrain.}

  \item{H2vec}{has weights for each RM relatednes matrix}

  \item{pheno}{A data frame with rownames set as sample IDs
  and a column containing phenotype data.}

  \item{phenoname}{The name of the column in pheno which
  contains phenotype data to test.}

  \item{Xcovamat}{Data frame of covariates with rownames()
  set to sample IDs.}

  \item{nfold}{Select the number of cross validation rounds
  to run. The value "LOOCV" will run one round of cross
  validation for each sample in your dataset.  The value
  "ncore" will set the test set size such that a single
  round runs on each core specified in the ncore option.
  Any numeric value will be set to the test size. Default
  runs 10 rounds of cross validation.}

  \item{ncore}{The number of cores available to distribute
  computaition across If a numeric value is supplied, that
  number of cores is registered. If the value "all" is
  supplied, all available cores are used.}

  \item{AUC}{Boolean value which sets whether or not a
  summary of the predictive predictive performance is
  generated.}
}
\value{
  A dataframe with three columns: sample ID, observed
  phenotype Ytest, and predicted phenotype Ypred
}
\description{
  This is a flexible cross validation routine which wraps
  the Omic Kriging calculation. The user can specify the
  size of the test set, all the way to "Leave One Out"
  cross validation. Additionally, all relevant parameters
  in the \code{\link{okriging}} function are exposed. This
  function uses the doMC package to distribute computation
  over multiple cores.
}
\keyword{prediction}

