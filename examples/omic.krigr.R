 ## Refactor the omic kriging using the GDS format

## dependencies
library(gdsfmt)
library(SNPRelate)
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
ncore <- 12

## load some functions for testing -- ordinarily you will simple load the package
source('R/gdsGeneRelMat.R')
source('R/dataInput.R')
source('R/computePCA.R')
source('R/krigrCrossValidation.R')

## load genetic data
load_gene_data(bedFile, bimFile, famFile, gdsFile)

## load phenotype data
pheno <- load_sample_data(phenoFile, main.pheno = pheno.name)

## calculate the genetic relatedness matrix
system.time(grm <- make_grm_gds(gdsFile = gdsFile, n.core = ncore))

## calculate principal components for use as covariates
system.time(pca <- make_PCs_irlba(grm, n.top = 2))


## n-fold parallel cross validation
result <- krigr_cross_validation(corlist = list(grm),
            pheno.df = pheno,
            pheno.name = pheno.name,
            Xcovmat = pca,
            H2vec = 1,
            ncore = 12,
            nfold = 10)




