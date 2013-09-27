#' Run n-fold cross validated kriging across n-cores.
#'
#' For now, this is a placeholder function. What is really needed is a function
#' which can run across multiple cores with an arbitrary number of test/train
#' groups, i.e. decouple the n-core and n-fold variables. At its core,
#' this function calls \code{\link{okriging}}.
#'
#' @param n.core The number of cores used and test sets the data will be split into.
#' @param corlist A list of correlation matrices. Colnames and rownames must be
#'   sample IDs, and must match betweenn matrices.
#' @param pheno.df A data frame of sample IDs and their corresponding phenotypes. The
#'   data frame must have rownames() set to sample IDs.
#' @param pheno.name The name of the column in the pheno.df which is the phenotype
#'   being predicted.
#' 
#' @return A data frame with a column for each: sample IDs, observed phenotypes,
#'   and predicted phenotypes.
#' 
#' @include R/okriging.R
#'
#' @keywords prediction, cross validation
#' @export
krigr_cross_validation <- function(n.cores, corlist, pheno.df, pheno.name) {
  ## TODO:: handling internal package references
  source('R/okriging.R')

  ## split into teating groups based on the number of cores available
  sample.ids <- pheno.df$IID
  n.samples <- length(sample.ids)
  groups <- 1:ncore
  rand.groups <- sample(groups, n.samples, replace=T)
  group.df <- data.frame(rand.groups, sample.ids)
  colnames(group.df) <- c("group.id", "sample.id")

  ## setup parallel environment
  registerDoMC(cores = ncore)

  ## running kriging routine on each core for each testing group
  p <- foreach(i=1:ncore, .combine=rbind) %dopar% {
    
    ## separate test/train for round i of the cross validation
    test.set = group.df$sample.id[group.df$group.id == i]
    train.set = group.df$sample.id[!(group.df$sample.id %in% test.set)]

    ## run kriging
    okriging(idtest=test.set,
      idtrain=train.set,
      corlist=corlist,
      H2vec=0.5, #?
             pheno=pheno.df,
             phenoname=pheno.name)
    }

  return(p)

}
