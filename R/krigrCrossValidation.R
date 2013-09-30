## n-core, n-fold cross validation routine for kriging
## TODO:: separate distinction between level of parallelism and fold level
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
