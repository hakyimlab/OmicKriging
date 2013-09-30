#' Multithreaded cross validation routine for Omic Kriging.
#'
#' This is a flexible cross validation routine which wraps the Omic Kriging
#' calculation. The user can specify the size of the test set, all the way to
#' "Leave One Out" cross validation. Additionally, all relevant  parameters in the
#' \code{\link{okriging}} function are exposed. This function uses the doMC
#' package to distribute computation over multiple cores.
#'
#' @param corlist A list of correlation matrices used in Kriging. rownames and colnames
#'   of cor should be IID list and include idtest and idtrain.
#' @param H2vec has weights for each RM relatednes matrix
#' @param pheno A data frame with rownames set as sample IDs and a column containing phenotype data.
#' @param phenoname The name of the column in pheno which contains phenotype data to test.
#' @param Xcovamat Data frame of covariates with rownames() set to sample IDs. 
#' @param nfold Select the number of cross validation rounds to run. The value "LOOCV"
#'   will run one round of cross validation for each sample in your dataset.
#'   The value "ncore" will set the test set size such that a single round
#'   runs on each core specified in the ncore option. Any numeric value
#'   will be set to the test size. Default runs 10 rounds of cross validation.
#' @param ncore The number of cores available to distribute computaition across
#'    If a numeric value is supplied, that number of cores is registered. If the
#'    value "all" is supplied, all available cores are used. 
#' @param AUC Boolean value which sets whether or not a summary of the predictive
#'    predictive performance is generated.
#'
#' @return  A dataframe with three columns: sample ID, observed phenotype Ytest, and predicted phenotype Ypred
#'
#' @keywords prediction
#'
#' @import doMC
#' @import ROCR
#' @export
krigr_cross_validation <- function(corlist, pheno.df, pheno.name, Xcovamat = NULL, H2vec, nfold = 10, ncore = "all", AUC = FALSE, ...) {
  ## TODO:: handling internal package references
  source('/nas40t0/keston/PS2/tempo/R/okriging.R')
  ## dependencies
  require(doMC)

  ## split into groups based on the number of cores available
  sample.ids <- pheno.df$IID
  n.samples <- length(sample.ids)

  ## detect cores  
  if(ncore == "all") {
    ncore <- detectCores()
    registerDoMC(cores = ncore)
    } else {
    registerDoMC(cores = ncore)
    }
  
  ## set n-fold  
  if(nfold == "LOOCV") {
    nfold <- n.samples
    } else if(is.numeric(nfold)) {
    nfold <- nfold
    } else if(nfold == "ncore") {
    nfold <- ncore
    } else {
    nfold <- 10
    }
   
  ## print core and fold numbers
  '%&%' <- function(a, b) paste(a, b, sep="")
  
  if(nfold == "LOOCV") {
    print('Set leave-one-out cross-validation...')
    } else {
    print('Set '%&% nfold %&%'x cross-validation...')
    }
  
  print('With '%&% ncore %&%' logical cores...') 
  
  ## create groups
  groups <- 1:nfold
  rand.groups <- sample(groups, n.samples, replace=T)
  group.df <- data.frame(rand.groups, sample.ids)
  colnames(group.df) <- c("group.id", "sample.id")

  print('Running OmicKriging...')
  
  ## running kriging routine on each core for each testing group
  time <- system.time(
  res <- foreach(i = 1:nfold, .combine = rbind) %dopar% {
    
    ## separate test/train for round i of the cross validation
    test.set <- group.df$sample.id[group.df$group.id == i]
    train.set <- group.df$sample.id[!(group.df$sample.id %in% test.set)]
    
    ## run kriging
    if(!is.null(Xcovamat)) {
      okriging(idtest = test.set, idtrain = train.set, corlist = corlist, 
      H2vec = H2vec, pheno = pheno.df, phenoname = pheno.name, Xcova = Xcovamat)
      } else {
      okriging(idtest = test.set, idtrain = train.set, corlist = corlist, 
      H2vec = H2vec, pheno = pheno.df, phenoname = pheno.name)

      }
    }
  )
 
  gc()

  ## summary
  if(AUC) {
    auc <- function(predtype,phenotype){
      require(ROCR)
      pred <- prediction(predtype,phenotype)
      perf <- performance(pred,"auc")
      aucval <- perf@y.values
      return(aucval)
      }
    print('Area under the ROC curve: '%&% auc(res$Ypred,res$Ytest))
    sum <- summary(lm(Ytest ~ Ypred, data = res))
    print(sum)
    } else {
    sum <- summary(lm(Ytest ~ Ypred, data = res))
    print(sum)
    }

  print('Finished OmicKriging in '%&% time[3] %&%' seconds')
  return(res)

}
