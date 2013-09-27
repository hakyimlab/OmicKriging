## n-core, n-fold cross validation routine for kriging
## TODO:: separate distinction between level of parallelism and fold level
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
