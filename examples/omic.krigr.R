### Refactor the omic kriging using the GDS format

## dependencies
library(gdsfmt)
library(SNPRelate)
library(Rcpp)
library(RcppEigen)
library(inline)
library(doMC)

## globals (can be factored out into high level config)
gdsFile <- "examples/test.gds"
bedFile <- "data/T1DCC.subset.bed"
bimFile <- "data/T1DCC.subset.bim"
famFile <-"data/T1DCC.subset.fam"
grmDataFile <- "data/T1DCC.subset.GRM.Rdata"
phenoFile <- "data/T1DCC.pheno"
pheno.name <- "PHENO"
ncore <- 4


## load genetic data
source('R/computeGeneRelMat.R')
source('R/data_loading.R')
load_gene_data(bedFile, bimFile, famFile, gdsFile)

## calculate the genetic relatedness matrix
grm <- make_grm(gdsFile = gdsFile)


## load phenotype data
pheno <- read.table(phenoFile, header=T)
## subset to individuals in GRM and match the sample order
sub.pheno.df <- subset( pheno.df, pheno.df$IID %in% colnames(grm) )
pheno <- sub.pheno.df
pheno[pheno.name] <- as.numeric(unlist(pheno[pheno.name]))
rownames(pheno) <- pheno$IID

## calculate principal components
pca <- make_PCs(gdsFile, n.core = ncore, n.top = 2)


## n-fold parallel cross validation
source('R/krigrCrossValidation.R')
result <- krigr_cross_validation(n.cores = ncore,
            corlist = list(grm),
            pheno.df = pheno,
            pheno.name = pheno.name)


## summarize results: prediction and performance functions from ROCR


