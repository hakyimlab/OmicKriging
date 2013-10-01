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
ReadGRMBin <- function(prefix, size = 4){
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
writeGRMBin <- function(X, n.snps = 0.0, prefix, size = 4) {

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

