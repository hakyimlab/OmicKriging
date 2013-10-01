#' Compute gene expression correlation matrix.
#'
#' This function computes a gene expression correlation matrix given a file of
#' transcript expression levels for each sample in the study. It returns a
#' correlation matrix with rownames and colnames as sample IDs.
#'
#' @param genefile Path to gene expression file.
#' @param outFilePrefix File path prefixes for outputting GCTA style binary
#'   correlation matrices.
#'
#' @return Returns a correlation matrix of (N-samples x N-samples), with
#'   rownames and colnames as sample IDs.
#'
#' @include R/grm_io.R
#'
#' @export
computeGX <- function(genefile, outFilePrefix, idfile=NULL) {
  
  source('R/grm_io.R')

  ## data input
  genedata = read.delim(genefile, sep="", as.is=T, header=T)
  if(!is.null(idfile)) {
    iddata = read.table(idfile,header=T,as.is=T)
    genedata = merge(iddata,genedata,by.x=c("FID","IID"),by.y=c("FID","IID"))
  }
  cordata.id = genedata[c("FID","IID")]
  genemat = as.matrix(genedata[,!(names(genedata) %in% c("FID","IID"))])
  ## impute missing data
  ## Standardize this matrix?
  if(sum(is.na(genemat))>0) print("Impute first") ## TODO: impute with the script
  cormat = cor(t(genemat))
  colnames(cormat) <- cordata.id[,2]
  rownames(cormat) <- cordata.id[,2]

  ## write to disk
  writeGRMBin(X = cormat, prefix = outFilePrefix, n.snps = length(genemat[1,])) 

  return(cormat)
}
