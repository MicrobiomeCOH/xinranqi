#Enable command line arguments
args <- commandArgs(TRUE)
theta <- as.numeric(args[1])
epsSD <- as.numeric(args[2])
numIndependentDiseased <- as.numeric(args[3])
numIndependentControl <- as.numeric(args[4])
numFamilies <- as.numeric(args[5])
numOfSelectedSnips <- as.numeric(args[6])
diseaseSnip <- as.numeric(args[7])
location <- as.character(args[8])
out_type <- as.character(args[9])

library(SKAT)
library(Matrix)
library(knockoff)
library(WGScan)
library(numDeriv)
library(glmnet)
library(SPAtest)
library(CompQuadForm)
library(SNPknock)
library(STAAR)
library(xtable)
library(ACAT)
#library(GhostKnockoff)
library(doMC)
library(irlba)
library(ff)
library(pROC)

#setwd('/scratch/users/xinranqi/Project1/continuous')
#setwd('/Users/xqi/Documents/Week2_01112022/family_simulation')
#source('/Users/xqi/Documents/Week6_02132022/KnockoffScreen_Nov5.R')
source('/oak/stanford/groups/zihuai/XinranQi/Project1/simulate10PedigreeRecombine0.R')
source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
library(KnockoffScreen)

data("SKAT.haplotypes")
SKATHaplotypes <- SKAT.haplotypes$Haplotype[ , SKAT.haplotypes$SNPInfo$FREQ1>0.001]
commonSnipPositionsInChroms <- SKAT.haplotypes$SNPInfo$CHROM_POS[SKAT.haplotypes$SNPInfo$FREQ1>0.001]
#commonSnipIndiceInPhasedFiles <- matrix(1:length(commonSnipPositionsInChroms), nrow=1)

phi <- matrix(c(1, 0, 0.5, 0.5, 0, 0, 0.25, 0.25, 0.25, 0.25,
                0, 1, 0.5, 0.5, 0, 0, 0.25, 0.25, 0.25, 0.25,
                0.5, 0.5, 1, 0.5, 0, 0, 0.5, 0.5, 0.25, 0.25,
                0.5, 0.5, 0.5, 1, 0, 0, 0.25, 0.25, 0.5, 0.5,
                0, 0, 0, 0, 1, 0, 0.5, 0.5, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0.5, 0.5,
                0.25, 0.25, 0.5, 0.25, 0.5, 0, 1, 0.5, 0.125, 0.125,
                0.25, 0.25, 0.5, 0.25, 0.5, 0, 0.5, 1, 0.125, 0.125,
                0.25, 0.25, 0.25, 0.5, 0, 0.5, 0.125, 0.125, 1, 0.5,
                0.25, 0.25, 0.25, 0.5, 0, 0.5, 0.125, 0.125, 0.5, 1), nrow=10, ncol=10, byrow=TRUE)

# Section 1, parameters for simulation
numMembersPerFamily <- rep(10, numFamilies)	#generate family sizes 10 members
heritageWeights <- 1 
H0DiseaseRates <- 0.2  #corresponding to the populations

#Section 2, under H1, define the disease rate associated with one snip
#rows corresponding to H1DiseaseModels: additive, dominate, recessive, multiplicative
#columns corresponding to the summation of the two snips of a pair of chromosomes
#00->0, 01->1, 10->1, 11->2
H1DiseaseRates <- matrix(c(0.1,0.15,0.20,  0.1,0.2,0.2,  0.1,0.1,0.2,  0.1,0.15,0.225),ncol=3,byrow=TRUE)

#under H1, which Snip is associated with disease?
#diseaseSnip <- sample(1:numOfSelectedSnips,1)  
#diseaseSnip <- 10

#one disease model chosen from the 4 models
#can be "HNull", "additive", "dominate", "recessive", "multiplicative"
thisDiseaseModel <- "recessive"
overallRecombinationRate <- 0 	
set.seed(123)

for (repIndx in 1:1000){
  pValues <- as.data.frame(matrix(NA, nrow=numOfSelectedSnips*2, ncol=3))
  pValues[ , 3] <- rep(c("STAAR", "ACAT"), each=numOfSelectedSnips)
  colnames(pValues) <- c("window", "single", "Type")
  
  if (repIndx %in% seq(1,1000,by=5)) {
    # select some snips to be included in the simulation
    selectedSnipIndiceInPhasedFiles <- sort(sample(ncol(SKATHaplotypes), size=numOfSelectedSnips, replace=FALSE))
    selectedSnipPositionsInChroms <- commonSnipPositionsInChroms[selectedSnipIndiceInPhasedFiles]
    
    # Section 5, read those selected snips from the phased files
    selectedSnipsFromPhasedFiles <- list()
    selectedSnipsFromPhasedFiles <- c(selectedSnipsFromPhasedFiles, list(SKATHaplotypes[,selectedSnipIndiceInPhasedFiles]))
    # prepare some global variables
    selectedSnipPositionRange <- range(commonSnipPositionsInChroms[selectedSnipIndiceInPhasedFiles])
    
    adjustedRecombinationRate <- overallRecombinationRate*(selectedSnipPositionRange[2]-selectedSnipPositionRange[1])/(max(commonSnipPositionsInChroms)-min(commonSnipPositionsInChroms))
    diseaseModels <- c("HNull","additive","dominate","recessive","multiplicative")
    isHNull <- thisDiseaseModel == diseaseModels[1]
    if (!isHNull) {
      for (i in 1:length(diseaseModels)) {
        if (thisDiseaseModel == diseaseModels[i]) {
          thisDiseaseRate <- H1DiseaseRates[i-1,]
        }
      }
    }
    tempPopMatrix <- matrix(rep(1:length(selectedSnipsFromPhasedFiles),numOfSelectedSnips),nrow=length(selectedSnipsFromPhasedFiles),numOfSelectedSnips)

    # generate related genotypes (in form of haplotypes) from 10 member pedigree
    chromosomes <- generateFamily()$a
    
    # convert to genotypes (coding 0,1,2)
    G <- matrix(0, nrow=dim(chromosomes)[1]/2, ncol=dim(chromosomes)[2])
    for (i in 1:(dim(chromosomes)[1]/2)){
      for (j in 1:dim(chromosomes)[2]){
        G[i,j] <- chromosomes[2*i-1,j] + chromosomes[2*i,j]
     }}
    
    for (j in 1:dim(chromosomes)[2]){
      if (!is.na(table(G[,j])[3])){ 
        if (table(G[,j])[1] < table(G[,j])[3]){
          #Switch 0 and 2 when number of 2 is larger than 0
          G[,j][G[,j]==0] <- 9
          G[,j][G[,j]==2] <- 0
          G[,j][G[,j]==9] <- 2
        }}
      if (!is.na(table(G[,j])[2])) {
        if (any(names(table(G[,j]))!=0) & is.na(table(G[,j])[3]) & names(table(G[,j])[2])==2){
          G[,j][G[,j]==2] <- 0
        }
      }
    }
    
    # random effect from kinship matrix & variance component parameter
    b1 <- matrix(NA, nrow=numFamilies, ncol=10)
    for (pedigreeSampleIndx in 1:numFamilies) {
      b1[pedigreeSampleIndx, ] <- mvrnorm(n=1, mu=rep(0, times=10), Sigma=theta*phi)
    }
    b1 <- as.vector(t(b1))
    
    # create kinship matrix for simulated observations
    kinshipMatrix <- matrix(0, nrow=length(b1), ncol=length(b1))
    for (pedigreeSampleIndx in 1:numFamilies) {
      kinshipMatrix[((pedigreeSampleIndx-1)*10+1):(pedigreeSampleIndx*10), ((pedigreeSampleIndx-1)*10+1):(pedigreeSampleIndx*10)] <- phi
    }
    colnames(kinshipMatrix) <- rownames(kinshipMatrix) <- 1:length(b1)
    # compressed sparse column format
    kinshipMatrix <- Matrix(kinshipMatrix,sparse=T)
    
    # restrict simulations to variants with minor allele counts >25 to ensure stable calculation of summary statistics
    MAC <- apply(G,2,sum)
    G <- G[,MAC>25]
    selectedSnipPositionsInChroms <- selectedSnipPositionsInChroms[MAC>25]
    
    # clustering
    cor.X <- cor(G)
    #cor.X <- sparse.cor(Matrix(G))$cor 
    Sigma.distance <- as.dist(1-abs(cor.X))
    fit <- hclust(Sigma.distance, method="single")
    corr_max <- 0.75
    clusters <- cutree(fit, h=1-corr_max)
    cluster.index <- match(unique(clusters), clusters)
    
    G <- G[ , cluster.index, drop=F]
    pos <- selectedSnipPositionsInChroms[cluster.index]

    Z <- as.matrix(rep(1,ncol(G)))
    window.size <- 2000
    weights <- as.matrix(dbeta(apply(G,2,mean)/2,1,25))#rep(1,ncol(SNP.set))#
    Z <- weights*Z
    #generate a matrix to specify the variants in each window
    window.matrix0 <- c()
    for(size in window.size){
      if (size==1) {next}
      pos.tag <- seq(min(pos), max(pos), by=size*1/2)
      pos.tag <- sapply(pos.tag, function(x) {pos[which.min(abs(x-pos))]})
      window.matrix0 <- cbind(window.matrix0, sapply(pos.tag, function(x) {as.numeric(pos>=x & pos<x+size)}))
    }
    
    #merge identical windows to get actual windows (start and end with variant position)
    window.string <- apply(window.matrix0, 2, function(x) {paste(as.character(x),collapse="")})
    window.matrix0 <- as.matrix(window.matrix0[ , match(unique(window.string),window.string)])
    window.matrix0 <- window.matrix0[ , apply(window.matrix0,2,sum)>1]
    #single variants will be added back later
    if(1 %in% window.size) {window.matrix0 <- as.matrix(window.matrix0[,apply(window.matrix0,2,sum)>1])} 
    #incorporate weights (MAF/annotations)
    window.matrix<-c()
    for(j in 1:ncol(Z)) {window.matrix <- cbind(window.matrix,Z[,j]*window.matrix0)}
    #calculate scan statistics for all windows
    window.summary <- t(apply(window.matrix0, 2, function(x) {c(min(pos[which(x==1)]),max(pos[which(x==1)]))}))
    # Generate the knockoff features: single variant
    window.matrix <- Matrix(window.matrix,sparse=T)
  }
  # simulate related phenotypes
  # under H0, all variants have no signal
  beta <- matrix(0, nrow=ncol(G), ncol=1)
  
  # whether the response variance is continuous 
  if (out_type == "C") {
    X1 <- rnorm(nrow(G),0,1)
    # plan 3: constant sum of theta
    epsSD <- ifelse(8>var(b1), sqrt(8-var(b1)), 0)
    mu <- X1+rnorm(nrow(G),0,epsSD)
    y <- as.matrix(mu+G%*%beta+b1)
  } else {
    prevalence <- 0.1
    b0_pheno <- log(prevalence/(1-prevalence)) # for binary trait
    X1 <- rnorm(nrow(G),0,1)
    mu <- X1+G%*%beta+b1
    mu <- mu-mean(mu)
    y <- as.matrix(rbinom(nrow(G),1, exp(mu+b0_pheno)/(1+exp(mu+b0_pheno))))
  }
  #y.pedigree <- matrix(y, nrow=nrow(G)/10, ncol=10, byrow=TRUE) 
  #corrPedigree <- cor(y.pedigree)
  #bmatrix(round(corrPedigree, digits=3))
  #pIndx <- as.vector(t(sapply(c(1,2,5,6), function(x){(1:400-1)*10+x })))
  #idx=50
  # correlation of phenotypes
  #x.pedigree <- matrix(G[pIndx, idx], nrow=nrow(G[pIndx, ])/4, ncol=4, byrow=TRUE) 
  #x.pedigree <- matrix(G[, idx], nrow=nrow(G)/10, ncol=10, byrow=TRUE) 
  #corrPedigree <- cor(x.pedigree)
  #bmatrix(round(corrPedigree, digits=3))
  # STAAR adjustment for kinship
  phenoSTAAR <- data.frame(Y=y,X1=X1,id=1:length(y))
  if (out_type=="C") {
    obj_nullmodel <- fit_null_glmmkin(Y~X1,data=phenoSTAAR,family=gaussian(link="identity"),id="id",kins=kinshipMatrix)
  } else {
    obj_nullmodel <- fit_null_glmmkin(Y~X1,data=phenoSTAAR,family=binomial(link="logit"),id="id",kins=kinshipMatrix)
  }
  # window based using STAAR
  pValuesSTAAR <- sapply(1:ncol(window.matrix), function(x) {STAAR.modified(genotype=G[ , window.matrix[ , x]!=0],obj_nullmodel=obj_nullmodel,rare_maf_cutoff=1,yResponse=y,out_type=out_type)$results_STAAR_O})
  # single variant using STAAR
  pSTAAR.add <- Indiv_Score_Test_Region(genotype=G, obj_nullmodel=obj_nullmodel, rare_maf_cutoff=1, rv_num_cutoff=1)[ , "pvalue"]

  # ACAT without adjustment for kinship
  # window based using ACAT
  obj <- NULL_Model(as.vector(y),Z=as.matrix(X1))
  pValuesACAT <- sapply(1:ncol(window.matrix), function(x) {ACAT_V(G[ , window.matrix[ , x]!=0],obj)})
  # single variant using ACAT 
  if (out_type == "C") {
    obj_nullmodel_identity <- fit_null_glm(Y~X1,data=phenoSTAAR,family=gaussian(link="identity"))
  } else {
    obj_nullmodel_identity <- fit_null_glm(Y~X1,data=phenoSTAAR,family=binomial(link="logit"))
  }
  pACAT.add <- Get.p.modifiedGLM(G=G, obj_nullmodel_glm=obj_nullmodel_identity)

  # fulfillment of p value matrix
  pValues[1:length(pValuesSTAAR), 1] <- pValuesSTAAR
  pValues[1:length(pSTAAR.add), 2] <- pSTAAR.add
  pValues[(1:length(pValuesACAT))+numOfSelectedSnips, 1] <- pValuesACAT 
  pValues[(1:length(pACAT.add))+numOfSelectedSnips, 2] <- pACAT.add 
  pValues <- as.data.frame(pValues[c(1:max(length(pValuesSTAAR), length(pSTAAR.add), length(pValuesACAT), length(pACAT.add)),
                                     1:max(length(pValuesSTAAR), length(pSTAAR.add), length(pValuesACAT), length(pACAT.add))+numOfSelectedSnips), ])
  # output
  write.table(pValues, paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/', location,'/ACATSTAAR_pValue_',out_type, '_',repIndx,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}
