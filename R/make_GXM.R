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
make_GXM <- function(expFile = NULL, gxmFilePrefix = NULL, idfile = NULL) {
  ## source file
  source('R/grm_IO.R')
  source('R/rcpp_CORMAT.R') 

  ## data input
  genedata <- read.delim(expFile, sep="", as.is=T, header=T)
  if(!is.null(idfile)) {
    iddata <- read.table(idfile,header=T,as.is=T)
    genedata <- merge(iddata,genedata,by.x=c("FID","IID"),by.y=c("FID","IID"))
  }
  cordata.id <- genedata[c("FID","IID")]
  genemat <- as.matrix(genedata[,!(names(genedata) %in% c("FID","IID"))])
  
  ## center and scale genemat
  genemat <- scale(genemat, center = TRUE, scale = TRUE)
  
  ## impute mean for missing values
  genemat[is.na(genemat)] <- 0.0
   
  ## compute cor mat
  cormat <- rcppcormat(t(genemat))
  
  ## give it row names
  colnames(cormat) <- cordata.id[,2]
  rownames(cormat) <- cordata.id[,2]

  ## write to disk
  if(!is.null(gxmFilePrefix)) {
    writeGRMBin(X = cormat, prefix = gxmFilePrefix, n.snps = length(genemat[1,])) 
  }
  return(cormat)
}
