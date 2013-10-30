#' Multithreaded cross validation routine for Omic Kriging.
#'
#' This is a flexible cross validation routine which wraps the Omic Kriging
#' calculation. The user can specify the size of the test set, all the way to
#' "Leave One Out" cross validation. Additionally, all relevant  parameters in the
#' \code{\link{okriging}} function are exposed. This function uses the doMC
#' package to distribute computation over multiple cores. If the phenotype is 
#' case/control, a ROCR AUC and GLM analysis is run and the results printed to screen.
#'
#' @param corlist A list of correlation matrices used in Kriging. rownames and colnames
#'   of cor should be IID list and include idtest and idtrain.
#' @param h2.vec has weights for each RM relatednes matrix
#' @param pheno.df A data frame with rownames set as sample IDs and a column containing phenotype data.
#' @param pheno.id The name of the column in pheno which contains phenotype data to test.
#' @param covar.mat Data frame of covariates with rownames() set to sample IDs. 
#' @param nfold Select the number of cross validation rounds to run. The value "LOOCV"
#'   will run one round of cross validation for each sample in your dataset.
#'   The value "ncore" will set the test set size such that a single round
#'   runs on each core specified in the ncore option. Any numeric value
#'   will be set to the test size. Default runs 10 rounds of cross validation.
#' @param ncore The number of cores available to distribute computaition across
#'    If a numeric value is supplied, that number of cores is registered. If the
#'    value "all" is supplied, all available cores are used. 
#' @param verbose Report rounds on cross validation on standard out. 
#'
#' @return  A dataframe with three columns: sample ID, observed phenotype Ytest, and predicted phenotype Ypred
#'
#' @keywords prediction, cross validation
#'
#' @include R/omic_KRIGR.R
#'
#' @import doMC
#' @import ROCR
#' @export
krigr_cross_validation <- function(cor.list, pheno.df, pheno.id = 1, h2.vec, covar.mat = NULL, nfold = 10, ncore = "all", verbose = FALSE, ...) {
  ## TODO:: handling internal package references
  source('R/omic_KRIGR.R')
  ## dependencies
  require(doMC)
  ## functions
  '%&%' <- function(a, b) paste(a, b, sep="")

  ## split into groups based on the number of cores available
  rownames(pheno.df) <- pheno.df$IID
  sample.ids <- pheno.df$IID
  n.samples <- length(sample.ids)
  cat("Detected", n.samples, "samples...", "\n")    
  flush.console()
  
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
  if(nfold == n.samples) {
    cat("Set leave-one-out cross-validation...", "\n")
    flush.console()
    } else {
    cat("Set", nfold %&% "x", "cross-validation...", "\n")
    flush.console()
    }
  
  cat("With", ncore, "logical core(s)...", "\n") 
  flush.console()

  ## create groups
  if(nfold == n.samples) {
    rand.groups <- 1:n.samples
    } else {
    rand.groups <- sample(1:nfold, n.samples, replace = T)
    }
  group.df <- data.frame(rand.groups, sample.ids)
  colnames(group.df) <- c("group.id", "sample.id") 
 
  cat("Running OmicKriging...", "\n")
  flush.console()
  
  ## running kriging routine on each core for each testing group
  n.par <- unique(rand.groups)
  time <- system.time(
    res <- foreach(i = 1:length(n.par), .combine = rbind) %dopar% {
      if (verbose) cat(n.par[i], "\n")
      flush.console()
    
      ## separate test/train for round i of the cross validation
      if(length(n.par) == n.samples){
        idtest <- as.character(sample.ids[i])
        } else {
        idtest <- as.character(group.df$sample.id[group.df$group.id == n.par[i]])
        }
      idtrain <- as.factor(group.df$sample.id[!(group.df$sample.id %in% idtest)])
      
      ## run kriging
      okriging(idtest, idtrain, 
            corlist <- cor.list, 
            H2vec <- h2.vec,
            pheno <- pheno.df,
            phenoname <- colnames(pheno.df)[pheno.id + 2],
            Xcova <- covar.mat
            )
      }
    )

  gc()

  ## summary
  if(length(unique(res$Ytest)) == 2) {
    auc <- function(predtype, phenotype){
      require(ROCR)
      pred <- prediction(predtype, phenotype)
      perf <- performance(pred, "auc")
      aucval <- perf@y.values
      return(aucval)
      }
    cat("Summary of binary phenotype...", "\n")
    cat("Area under the ROC curve:", auc(res$Ypred,res$Ytest) %&% "...", "\n")
    convertpheno <- res$Ytest
    uniconvertpheno <- unique(convertpheno)
    convertpheno[convertpheno == uniconvertpheno[2]] = 0
    convertpheno[convertpheno == uniconvertpheno[1]] = 1
    sum <- summary(glm(convertpheno ~ res$Ypred, family = binomial))
    print(sum)
    } else {
    sum <- summary(lm(Ytest ~ Ypred, data = res))
    print(sum)
    }

  cat("Finished OmicKriging in", time[3], "seconds", "\n")
  return(res)

}
