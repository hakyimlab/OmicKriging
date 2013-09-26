### Refactor the omic kriging using the GDS format

## all dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)
library(doMC)

## functions


## globals to be factored out into high level config
gdsFile <- "wtcct1d.gds"
bedFile <- "data/T1DCC.bed"
bimFile <- "data/T1DCC.bim"
famFile <-"data/T1DCC.fam"
ncore <- 12

## INPUTS 
## convert plink file to GDS
snpgdsBED2GDS(bedFile, famFile, bimFile, gdsFile)
## open GDS file and pull genotype matrix
genofile <- openfn.gds(gdsFile)
X <- read.gdsn(index.gdsn(genofile, "genotype"))
Y <- read.gdsn(index.gdsn(genofile, "sample.annot"))$phenotype

## TODO:: factor phenotype handling into GDS structure
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
pheno.df <- dataframe(test=Y)
rownames(pheno.df) <- sample.id

## compute GRM matrix (or read in?)
source('src/rcppcormat.r')
## center genotype matrix
X.center <- sweep(X, 2, colMeans(X), "-")
grm <- rcppcormat(t(X.center))

## PARALLEL KRIGING
source('src/okriging.R')
corlist <- list()
corlist[[1]] <- grm


## setup the cores
registerDoMC(cores = ncore)

## parallel cross validation (n-fold)
idlength <- length(sample.id)
g <- 1:ncore
groupid <- sample(g,idlength,replace=T)
newiddata <- data.frame(groupid, sample.id)
colnames(newiddata) <- c("group.id","sample.id")

## run kriging in parallel
p <- foreach(i=1:ncore, .combine=rbind) %dopar% {
  
  test.set = newiddata$sample.id[newiddata$group.id == i]
  train.set = newiddata$sample.id[!(newiddata$sample.id %in% test.set)]
  
  okriging(idtest=test.set,
           idtrain=train.set,
           corlist=corlist,
           H2vec=0.5, #?
           pheno=pheno.df,
           phenoname="test")
}





