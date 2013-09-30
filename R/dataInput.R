#' Load genetic data from PLINK binary files.
#'
#' This function takes file paths for PLINK binary files and loads them into
#' a Genomic Data Structure (GDS) stored on disk. This is a wrapper function
#' aroud utilities provided by the gdsfmt and SNPRelate packages. The GDS is
#' used in downstream analyses such as: computing the genetic related matrix
#' (\code{\link{make_grm}}) and computing principal components (\code{\link{make_PCs}})
#'
#' @param bedFile File path to PLINK .bed file.
#' @param bimFile File path to PLINK .bim file.
#' @param famFile File path to PLINK .fam file.
#' @param gdsFile File path to store GDS for other analyses.
#'
#' @return None. The GDS file is stored on disk.
#'
#' @import gdsfmt
#' @import SNPRelate
#'
#' @keywords input
#'
#' @references library(gdsfmt), library(SNPRelate)
#'
#' @export
load_gene_data <- function(bedFile, bimFile, famFile, gdsFile) {
  require(gdsfmt)
  require(SNPRelate)
  snpgdsBED2GDS(bedFile, famFile, bimFile, gdsFile)
}

#' Loads sample phenotype and covariate data into data frame.
#'
#' This function loads a file into a data frame. This file should contain one
#' row per sample in your study, and one column for each covariate and
#' phenotype of interest. Additionally, it requires a header with "IID" for
#' the column of sample IDs, and a unique name for each phenotype and covariate.
#'
#' @param phenoFile File path to the phenotype/covariate file.
#' @param main.pheno Column name of the main phenotype of interest.
#'
#' @return A data frame with dimensions (# of samples) x (# of phenotypes/covar)
#'
#' @keywords input
#'
#' @export
load_sample_data <- function(phenoFile, main.pheno) {

  ## TODO:: make columns numeric when they can reasonably be converted into numeric

  # load phenotype data
  pheno <- read.table(phenoFile, header=T)
  pheno[main.pheno] <- as.numeric(unlist(pheno[main.pheno]))
  rownames(pheno) <- pheno$IID

  return(pheno)
}

