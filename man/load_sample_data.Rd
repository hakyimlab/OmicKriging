\name{load_sample_data}
\alias{load_sample_data}
\title{Loads sample phenotype and covariate data into data frame.}
\usage{
  load_sample_data(phenoFile, main.pheno)
}
\arguments{
  \item{phenoFile}{File path to the phenotype/covariate
  file.}

  \item{main.pheno}{Column name of the main phenotype of
  interest.}
}
\value{
  A data frame with dimensions (# of samples) x (# of
  phenotypes/covar)
}
\description{
  This function loads a file into a data frame. This file
  should contain one row per sample in your study, and one
  column for each covariate and phenotype of interest.
  Additionally, it requires a header with "IID" for the
  column of sample IDs, and a unique name for each
  phenotype and covariate.
}
\keyword{input}

