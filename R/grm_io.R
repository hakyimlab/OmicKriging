#' Read the GRM binary file.
#'
#' Function provided partially by GCTA maintainers for accessing their binary
#' GRM format introduced recently. The GRM is stored as a vector of numerics
#' which correspond to the lower triangular elements. We simply read these, pull
#' the diagonal elements, and inflate this into a full symmetric matrix. We add
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
  
  ## clean up file connections
  close(BinFile)
  close(NFile)

  return(X)
}

#' Write GRM binary files
#'
#'
#'
#'
#'
#'
writeGRMBin <- function {



}

