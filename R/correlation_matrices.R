#' Compute genetic correlation matrix from PLINK binary files.
#'
#' This is a convenience function which produces a centered genetic correlation
#' matrix from SNPs loaded into a Genomic Data Structure (GDS) file. The resulting matrix can be used 
#' with the okriging function. The GRM can be saved to disk as a
#' R object for fast loading downstream. The genotype  matrix is z-score normalized (i.e.
#' column means are centered and column variance is divided out to unit variance)
#' prior to calculating the correlation matrix.
#'
#' @param gdsFile File holding the GDS from which to pull the raw genotype matrix.
#' @param grmFilePrefix File to store the resulting GRM on disk as an R object.
#' @param snpList A vector of SNP IDs to subset the GRM on.
#' @param sampleList A vector of sample IDs to subset the GRM on.
#'
#' @return A genetic correlation matrix with colnames and rownames set to sample IDs.
#'   Each entry in the matrix is of type 'double'.
#'
#' @import gdsfmt
#' @import SNPRelate
#'
#' @keywords input, GRM
#' @export
make_GRM <- function(gdsFile = NULL, grmFilePrefix = NULL, snpList = NULL, sampleList = NULL) {

  genofile <- openfn.gds(gdsFile)
  ## pull an integer dosage matrix from the GDS. Rows are samples, columns are SNPs, and missing values are int 3.
  X <- snpgdsGetGeno(gdsobj = genofile, sample.id = sampleList, snp.id = snpList, verbose = FALSE)
  ## set missing values (int 3) to properly missing
  X[X == 3] <- NA
  ## z-normalize matrix (sweep out column means, and divide out column matrices)
  X <- scale(X, center = TRUE, scale = TRUE)
  ## set missing values to new column mean, i.e. 0.0
  X[is.na(X)] <- 0.0
  grm <- rcppcormat(t(X))
  
  ## pull sample IDs unless a sample list is specified
  if(!is.null(sampleList)) {
    sample.ids <- sampleList
  } else {
    sample.ids <- read.gdsn(index.gdsn(genofile, "sample.id"))
  }

  ## annotate columns and rows with sample IDs
  colnames(grm) <- sample.ids
  rownames(grm) <- sample.ids

  ## write out the GRM if a file is specified
  if(!is.null(grmFilePrefix)) {
    write_GRMBin(X = grm, prefix = grmFilePrefix)
  }

  return(grm)
}

#' Compute gene expression correlation matrix.
#'
#' This function computes a gene expression correlation matrix given a file of
#' transcript expression levels for each sample in the study. It returns a
#' correlation matrix with rownames and colnames as sample IDs.
#'
#' @param expFile Path to gene expression file.
#' @param gxmFilePrefix File path prefixes for outputting GCTA style binary
#'   correlation matrices.
#' @param idfile Path to file containing family IDs and sample IDs with header FID and IID. 
#'
#' @return Returns a correlation matrix of (N-samples x N-samples), with
#'   rownames and colnames as sample IDs.
#'
#' @export
make_GXM <- function(expFile = NULL, gxmFilePrefix = NULL, idfile = NULL) {

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
    write_GRMBin(X = cormat, prefix = gxmFilePrefix, n.snps = length(genemat[1,])) 
  }
  return(cormat)
}

#' Run Principal Component Analysis (PCA) using the Genomic Data Structure (GDS).
#'
#' An efficient method for computing Principal Components using the Genomic Data
#' Structure (GDS). This is a convenience wrapper for functions from the 
#' SNPRelate package.
#'
#' @param gdsFile A Genomic Data Structure file describing your study.
#' @param n.top Number of top principal components to return. Defaults to
#'   returning all components (i.e. # of samples).
#' @param n.core Distrubute computation across N cores.
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top).
#'
#' @keywords covariate, PCA, GRM
#'
#' @import SNPRelate
#'
#' @export
make_PCs_gds <- function(gdsFile, n.core, n.top = 0) {

  gds <- openfn.gds(gdsFile)
  pca <- snpgdsPCA(gds, num.thread = n.core)
  rownames(pca$eigenvect) <- read.gdsn(index.gdsn(gds, "sample.id"))
  return( pca$eigenvect )
}

#' Run Principal Component Analysis (PCA) using base R svd() function.
#'
#' A simple wrapper around the base R svd() function which returns the top N
#' eigenvectors of a matrix. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#' @export
make_PCs_svd <- function(X, n.top = 2) {
  res <- La.svd(X, nu = n.top)
  rownames(res) <- rownames(X)
  return(res["u"])
}

#' Run Principal Component Analysis (PCA) using the irlba package.
#'
#' A simple wrapper around the irlba() function which computes a partial SVD
#' efficiently. This function's run time depends on the number of eigenvectors
#' requested but scales well. Use this function to generate covariates for use
#' with the \code{\link{okriging}} or \code{\link{krigr_cross_validation}}
#' functions.
#'
#' @param X A correlation matrix.
#' @param n.top Number of top principal compenents to return 
#'
#' @return A matrix of Principal Components of dimension (# of samples) x
#'   (n.top). As expected, eigenvectors are ordered by eigenvalue. Rownames
#'   are given as sample IDs.
#'
#' @keywords covariate, PCA, GRM
#'
#' @references library(irlba)
#'
#' @import irlba
#' @export
make_PCs_irlba <- function(X, n.top = 2) {
  
  res <- irlba(X, nu = n.top)
  rownames(res$u) <- rownames(X)
  return(res$u)
}

#' Compute a correlation matrix.
#'
#' This function computes a correlation matrix using a cross product. It is
#' implemented in C++ for performance. This function is compiled and called
#' from R using the Rcpp package.
#'
#' @param snpmat A centered SNP matrix containing fields of type 'double'.
#'
#' @return An n x n correlation matrix.
#'
#' @import Rcpp
#' @import inline
#'
#' @keywords c++, performance, correlation matrix
rcppcormat <- function(snpmat){
    ## rcpp
    crossprodCpp <- '
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::Lower;

    // set constants

    const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
    const int n(A.cols());

    // create n x n matrix

    MatrixXd AtA(MatrixXd(n, n)

    // zero out the matrix

    .setZero()

    // treat what follows as symmetric matrix, use lower

    .selfadjointView<Lower>()

    // sum of B + AA

    .rankUpdate(A.adjoint()));

    // return

    return wrap(AtA);
    '
  
    ## compile the cross product function
    cpcpp <- cxxfunction(signature(AA="matrix"), crossprodCpp,
    plugin="RcppEigen", verbose=FALSE)

    cormat <- cpcpp(snpmat)

    ## post work
    cormatdiag <- diag(cormat)
    cormat <- sweep(cormat, 1, cormatdiag, "/")
    cormat <- sweep(cormat, 2, cormatdiag, "/")
    return(cormat)

}

#' Read the GRM binary file.
#'
#' Function provided by GCTA maintainers (modified slightly) for accessing their
#' recently introduced binary GRM format. The GRM is stored as a vector of numerics
#' which correspond to the lower triangular elements including the diagonal. We simply read these, pull
#' the diagonal elements, and inflate them into a full symmetric matrix. We add
#' sample IDs to colnames and rownames for compatibility with other Kriging 
#' functions.
#'
#' @param prefix The file path prefix to GRM binary files (e.g., test.grm.bin, test.grm.N.bin, test.grm.id.)
#' @param size The length (in bytes) of each value in the raw GRM vector. Default is 4, and matches GRM writen by GCTA 1.11.
#'
#' @return GRM of dim (N.samples x N.samples) with rownames and colnames as sample ID.
#'
#' @references http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
#'
#' @export
read_GRMBin <- function(prefix, size = 4){
  sum_i <- function(i){
    return(sum(1:i))
  }

  ## open file connections and read in data
  BinFileName <- paste(prefix,".grm.bin",sep="")
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  id <- read.table(IDFileName)
  n <- dim(id)[1]
  BinFile <- file(BinFileName, "rb")
  grm <- readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile <- file(NFileName, "rb")

  ## read in the number of SNPs used to calculate the GRM (does not appear to work)
  N <- readBin(NFile, n=1, what=numeric(0), size=size)
  i <- sapply(1:n, sum_i)
  
  ## clean up file connections
  close(BinFile)
  close(NFile)

  ## pull apart diagonal and lower triagular elements
  diag.elem <- grm[i]
  off.diag.elem <- grm[-i]

  ## create the full symmetric correlation matrix
  X <- diag(diag.elem)
  X[ lower.tri(X, diag = FALSE) ] <- off.diag.elem
  X <- X + t(X) - diag(diag(X)) 

  ## add sample IDs to rownames and colnames
  rownames(X) <- id$V2
  colnames(X) <- id$V2

  return(X)
}

#' Write GRM binary files.
#'
#' Function to write a binary GRM format recently introduced by GCTA. It takes
#' a correlation matrix as used by other Kriging functions, and writes three
#' files: binary file for storing the diagonal + lower triangular elements, a
#' text file for sample IDs, and a binary file storing the number of SNPs used
#' in the correlation matrix calculation.
#'
#' @param X Correlation matrix with rownames and colnames as sample IDs.
#' @param prefix Base file path and names for the three output files.
#' @param n.snps Number of SNPs used in correlation matrix calculation. Default is 0.0.
#' @param size Number of bytes to write for each value. Default is 4
#'
#' @return None. Though side effects are writing three files as described above.
#'
#' @references http://www.complextraitgenomics.com/software/gcta/estimate_grm.html
#'
#' @export
write_GRMBin <- function(X, n.snps = 0.0, prefix, size = 4) {

  sum_i <- function(i){
    return(sum(1:i))
  }

  ## file connections
  NFileName <- paste(prefix,".grm.N.bin",sep="")
  IDFileName <- paste(prefix,".grm.id",sep="")
  BinFileName <- paste(prefix,".grm.bin",sep="")

  ## pull sample ids and dimension of GRM
  id <- rownames(X)
  n <- length(id)

  ## pull diagonal elements
  diag.elem <- diag(X)

  ## pull lower triangular elements
  off.diag.elem <- X[lower.tri(X, diag=FALSE)]

  ## collapse GRM into vector
  i <- sapply(1:n, sum_i)
  collapsed.grm <- vector(mode="numeric", length = n*(n+1)/2)
  collapsed.grm[i] <- diag.elem
  collapsed.grm[-i] <- off.diag.elem

  ## write binary files
  BinFile <- file(BinFileName, "wb")
  NFile <- file(NFileName, "wb")
  writeBin(con = BinFile, collapsed.grm, size = size)
  writeBin(con = NFile, n.snps, size = size )
  close(BinFile)
  close(NFile)

  ## write sample ID file -- we are dropping sample family IDs here
  write.table(cbind(id, id), file = IDFileName)

}

