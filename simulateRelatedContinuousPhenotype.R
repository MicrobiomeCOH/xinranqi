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
h2 <- as.numeric(args[9])
N.causal <- as.numeric(args[10])
out_type <- as.character(args[11])

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
#diseaseSnip <- 10

#one disease model chosen from the 4 models
#can be "HNull", "additive", "dominate", "recessive", "multiplicative"
thisDiseaseModel <- "recessive"
overallRecombinationRate <- 0	

singleResult.detail <- c()
set.seed(123)

for (replicateIndx in 1:1000){
  if (replicateIndx %in% seq(1,1000,by=5)) {
    #Section 4, select some snips to be included in the simulation
    selectedSnipIndiceInPhasedFiles <- sort(sample(ncol(SKATHaplotypes), size=numOfSelectedSnips, replace=FALSE))
    selectedSnipPositionsInChroms <- commonSnipPositionsInChroms[selectedSnipIndiceInPhasedFiles]
    
    #Section 5, read those selected snips from the phased files
    selectedSnipsFromPhasedFiles <- list()
    selectedSnipsFromPhasedFiles <- c(selectedSnipsFromPhasedFiles, list(SKATHaplotypes[,selectedSnipIndiceInPhasedFiles]))
    #prepare some global variables
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
    weights <- as.matrix(dbeta(apply(G,2,mean)/2,1,25)) #rep(1,ncol(SNP.set))#
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

    # knockoff counterparts
    G_k <- create.MK(G,pos,M=5,maxN.neighbor=Inf,maxBP.neighbor=0.1*10^6)
    G_k <- if(is.list(G_k)) {G_k} else {list(G_k[1,,], G_k[2,,], G_k[3,,], G_k[4,,], G_k[5,,])}
    #G_k <- create.MK(G,pos,M=5,maxN.neighbor=Inf,maxBP.neighbor=0.1*10^6)
    G_k2 <- create.second_order(G,method="sdp",shrink=T)
    # HMM model
    G <- apply(G,2,as.integer)
    numOfVariant <- ncol(G) # Number of variables in the model
    numOfK <- 5 # Number of possible states for each variable
    M <- 3
    # Marginal distribution for the first variable
    pInit <- rep(1/numOfK,numOfK)
    # Create numOfVariant-1 transition matrices
    Q <- array(stats::runif((numOfVariant-1)*numOfK*numOfK),c(numOfVariant-1,numOfK,numOfK))
    for(j in 1:(numOfVariant-1)) {
      #Q[j,,] = Q[j,,] + diag(rep(1,numOfK))
      Q[j,,] <- Q[j,,] / rowSums(Q[j,,])
    }
    pEmit <- array(stats::runif(numOfVariant*M*numOfK),c(numOfVariant,M,numOfK))
    for(j in 1:numOfVariant) { 
      pEmit[j,,] <- pEmit[j,,] / rowSums(pEmit[j,,]) 
    }
    G_k3 <- knockoffHMM(G,pInit,Q,pEmit)

    #summary stats -- always use shrinkage
    cor.G <- matrix(as.numeric(corpcor::cor.shrink(G,verbose=F)), nrow=ncol(G))
    fit.prelim.M <- ES.prelim(cor.G,M=5) #fit preliminary model
    # fit null model
    fit.prelim.Ghost <- GhostKnockoff.prelim(cor.G,M=5,method='sdp',corr_max=0.75)
    
    #multiple second order knockoffs
    mu.G <- colMeans(G)
    sd.G <- apply(G,2,sd)
    scale.G <- t((t(G)-mu.G)/sd.G)
    E <- matrix(rnorm(nrow(scale.G)*5*ncol(scale.G)),nrow(scale.G),5*ncol(scale.G))
    scale.G_k2.M <- t(apply(scale.G%*%t(fit.prelim.Ghost$P.each),1,rep,times=5))+E%*%t(fit.prelim.Ghost$V.left)
    temp.G_k2.M <- t(rep(sd.G,5)*t(scale.G_k2.M)+rep(mu.G,5))
    
    # design matrices
    X <- as.matrix(G%*%window.matrix)
    X_k <- lapply(G_k,function(x) as.matrix(x)%*%window.matrix)
    X_k2 <- as.matrix(G_k2%*%window.matrix)
    X_k3 <- as.matrix(G_k3%*%window.matrix)
  }
  
  # simulate related phenotypes
  #signal.region<-10000
  #start <- sample(pos,1); end <- start+signal.region
  #signal.index <- which(pos>=start & pos<=end)
  #causal.index <- (1:ncol(G)) %in% sample(signal.index,ceiling(length(signal.index)*P.causal))
  causal.index <- (1:ncol(G)) %in% sample.int(ncol(G), size=N.causal, replace=FALSE)
  causal.index[sample(which(causal.index>0), N.causal/2, replace=FALSE)] <- -1
  # regression coefficients
  beta <- sqrt(h2/N.causal/apply(G,2,var))*causal.index
  
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
  
  # apply other alternative methods
  result.prelim <- KS.prelim(y,X=X1,out_type=out_type)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  
  #permute the residuals
  B<-1000
  permute.index<-sapply(1:B,function(x)sample(1:length(result.prelim$Y)))
  temp.Y.res<-Y.res[as.vector(permute.index)]
  re.Y.res<-matrix(temp.Y.res,length(result.prelim$Y),B)
  #re.Y.res<-result.prelim$re.Y.res
  score<-t(G)%*%Y.res;re.score<-t(t(G)%*%re.Y.res)
  #score_k<-apply(G_k,1,function(x) t(x)%*%Y.res);re.score_k<-lapply(1:dim(G_k)[1],function(s) t(t(G_k[s,,])%*%re.Y.res))
  score_k<-as.data.frame(lapply(G_k,function(x) t(as.matrix(x))%*%Y.res))
  re.score_k<-lapply(G_k,function(s) t(t(as.matrix(s))%*%re.Y.res))
  
  #all Cauchy test
  MAF<-apply(G,2,mean)/2
  MAC<-apply(G,2,sum)
  weight.beta<-dbeta(apply(G,2,mean)/2,1,25)
  weight.matrix<-cbind(MAC<10,(MAF<0.01&MAC>10)*weight.beta,(MAF>=0.01)*weight.beta)
  colnames(weight.matrix)<-c('MAC<5','MAF<0.01&MAC>5&Beta','MAF>=0.01Beta')
  #KS.fit<-KS.test(G,G_k,result.prelim,window.matrix=window.matrix0,weight.matrix=weight.matrix)
  p.add<-as.matrix(Get.p(G,result.prelim))
  p.add_k<-as.matrix(as.data.frame(lapply(G_k,function(x) Get.p(as.matrix(x),result.prelim=result.prelim))))
  # single variable Cauchy test
  #p.single<-KS.fit$p.single
  #p.single_k<-KS.fit$p.single_k
  #p.add<-as.matrix(p.single); p.add_k<-p.single_k
  KS.fit<-KS.test(temp.G=G,temp.G_k=G_k,p.rare=p.add,p.rare_k=p.add_k,result.prelim=result.prelim,window.matrix=window.matrix0,weight.matrix=weight.matrix)
  
  # SummaryStat knockoff - score test
  # knockoff screen - SCIT
  n<-length(Y.res)
  Zscore_0<-qnorm(p.add/2,lower.tail=F)*sign(cov(as.matrix(G),Y.res)) #convert p-values to z-scores
  ES.stat<-ES.fit(Zscore_0,n,gamma=1,fit.prelim=fit.prelim.M) #fit summary stat based knockoff
  MK.stat<-MK.statistic(ES.stat$T_0,ES.stat$T_k) #compute knockoff statistics
  W_summaryStat_score<-MK.q.byStat(MK.stat[,'kappa'],MK.stat[,'tau'],M=5,Rej.Bound=10000) #calculate q value
  
  # Ghostknockoff score test
  GK.stat<-GhostKnockoff.fit(Zscore_0=Zscore_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  #GK.stat <- GhostKnockoff.fit(Zscore_0,n,fit.prelim=fit.prelim.Ghost,gamma=1,weight.study=NULL)
  GK.filter <- MK.statistic(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
  #GK.filter <- GhostKnockoff.filter(GK.stat$T_0,GK.stat$T_k)
  #W_ghostScore <- rep(1,nrow(GK.stat$T_0))
  W_ghostScore <- MK.q.byStat(GK.filter[,'kappa'],GK.filter[,'tau'],5,Rej.Bound=10000)
  #W_ghostScore <- GK.filter$q
  
  p.add_k2 <- as.matrix(Get.p(G_k2,result.prelim=result.prelim))
  p.add_k3 <- as.matrix(Get.p(G_k3,result.prelim=result.prelim))
  
  # STAAR adjustment for kinship
  phenoSTAAR <- data.frame(Y=y,X1=X1,id=1:length(y))
  if (out_type=="C") {
    obj_nullmodel <- fit_null_glmmkin(Y~X1,data=phenoSTAAR,family=gaussian(link="identity"),id="id",kins=kinshipMatrix)
  } else {
    obj_nullmodel <- fit_null_glmmkin(Y~X1,data=phenoSTAAR,family=binomial(link="logit"),id="id",kins=kinshipMatrix)
  }
  # single variant using STAAR
  pSTAAR.add <- matrix(Indiv_Score_Test_Region(genotype=G, obj_nullmodel=obj_nullmodel, rare_maf_cutoff=1, rv_num_cutoff=1)[ , "pvalue"], ncol=1)
  pSTAAR.add_k <- matrix(NA, nrow=ncol(G), ncol=length(G_k))
  for (repIndx in 1:length(G_k)) {
    pSTAAR.add_k[ , repIndx] <- Indiv_Score_Test_Region(genotype=as.matrix(G_k[[repIndx]]), obj_nullmodel=obj_nullmodel, rare_maf_cutoff=1, rv_num_cutoff=1)[ , "pvalue"]
  }
  
  # multiple second order, STAAR
  pSecondOrderSTAAR.add_k <- matrix(NA, nrow=length(pSTAAR.add), ncol=5)
  for (secdIndx in 1:ncol(pSecondOrderSTAAR.add_k)) {
    secondG_k2 <- as.matrix(temp.G_k2.M)[ , ((secdIndx-1)*ncol(G)+1):(secdIndx*ncol(G))]
    pSecondOrderSTAAR.add_k[ , secdIndx] <- Indiv_Score_Test_Region(genotype=secondG_k2, obj_nullmodel=obj_nullmodel, rare_maf_cutoff=1, rv_num_cutoff=1)[ , "pvalue"]
  }

  # SummaryStat knockoff - mixed model score test
  # knockoff screen - SCIT
  n<-length(obj_nullmodel$residuals)
  ZscoreAdd_0<-qnorm(pSTAAR.add/2,lower.tail=F)*sign(cov(as.matrix(G),obj_nullmodel$residuals)) #convert p-values to z-scores
  ESSTAAR.stat<-ES.fit(ZscoreAdd_0,n,gamma=1,fit.prelim=fit.prelim.M) #fit summary stat based knockoff
  MKSTAAR.stat<-MK.statistic(ESSTAAR.stat$T_0,ESSTAAR.stat$T_k) #compute knockoff statistics
  W_summaryStat_STAAR<-MK.q.byStat(MKSTAAR.stat[,'kappa'],MKSTAAR.stat[,'tau'],M=5,Rej.Bound=10000) #calculate q value
  
  # Ghostknockoff mixed model score test
  GK.stat.STAAR <- GhostKnockoff.fit(Zscore_0=ZscoreAdd_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  #GK.stat.STAAR <- GhostKnockoff.fit(ZscoreAdd_0,n,fit.prelim=fit.prelim.Ghost,gamma=1,weight.study=NULL)
  GK.filter.STAAR <- MK.statistic(GK.stat.STAAR$T_0[[1]],GK.stat.STAAR$T_k[[1]])
  #GK.filter.STAAR <- GhostKnockoff.filter(GK.stat.STAAR$T_0,GK.stat.STAAR$T_k)
  #W_ghostScore.STAAR <- rep(1,nrow(GK.stat.STAAR$T_0))
  W_ghostScore.STAAR <- MK.q.byStat(GK.filter.STAAR[,'kappa'],GK.filter.STAAR[,'tau'],5,Rej.Bound=10000)
  #W_ghostScore.STAAR <- GK.filter.STAAR$q

  # proposed method: adjustment for sample relatedness while generating knockoffs
  # variance component parameter estimates from STAAR::fit_null_glmmkin
  thetaHat <- as.numeric(ifelse(is.list(obj_nullmodel$theta), obj_nullmodel$theta$kins1, obj_nullmodel$theta["kins1"]))
  #thetaHat <- as.numeric(obj_nullmodel$theta$kins1)
  # The conditional variance of reponse in Gaussian family is a constant: \phi (dispersion parameter), which is the alternative notation for the error variance: \sigma^2(e).
  dispersionHat <- as.numeric(ifelse(is.list(obj_nullmodel$theta), obj_nullmodel$theta$residuals, obj_nullmodel$theta["dispersion"]))
  #dispersionHat <- as.numeric(obj_nullmodel$theta$residuals)
  # check sum of variance == ?!& thetaHat diff theta 
  etaHat <- etaEstimate(N=n, ThetaHat=thetaHat, Var.e.Hat=dispersionHat, KinshipMatrix=kinshipMatrix)
  fit.prelimAdjusted.M<-ES.prelim.modified(cor.G,M=5,EtaHat=etaHat) #fit proposed adjustment for sample relatedness model
  
  ESSTAARAdjust.stat<-ES.fit(ZscoreAdd_0,n,gamma=1,fit.prelim=fit.prelimAdjusted.M) 
  MKSTAARAdjust.stat<-MK.statistic(ESSTAARAdjust.stat$T_0,ESSTAARAdjust.stat$T_k) 
  W_summaryStat_STAARAdjust<-MK.q.byStat(MKSTAARAdjust.stat[,'kappa'],MKSTAARAdjust.stat[,'tau'],M=5,Rej.Bound=10000)
  
  # ACAT without adjustment for kinship
  # window based using ACAT
  obj <- NULL_Model(as.vector(y),Z=as.matrix(X1))
  # single variant using ACAT 
  if (out_type == "C") {
    obj_nullmodel_identity <- fit_null_glm(Y~X1,data=phenoSTAAR,family=gaussian(link="identity"))
  } else {
    obj_nullmodel_identity <- fit_null_glm(Y~X1,data=phenoSTAAR,family=binomial(link="logit"))
  }
  pACAT.add <- matrix(Get.p.modifiedGLM(G=G, obj_nullmodel_glm=obj_nullmodel_identity), ncol=1)
  pACAT.add_k <- matrix(NA, nrow=ncol(G), ncol=length(G_k))
  for (repIndx in 1:length(G_k)) {
    pACAT.add_k[ , repIndx] <- Get.p.modifiedGLM(G=as.matrix(G_k[[repIndx]]), obj_nullmodel_glm=obj_nullmodel_identity)
  }

  pACAT.add_secondOrder <- matrix(Get.p.modifiedGLM(G=G_k2, obj_nullmodel_glm=obj_nullmodel_identity), ncol=1)
  # multiple second order, GLM single variant
  pACAT.add_secondOrder.M <- matrix(NA, nrow=length(pACAT.add_secondOrder), ncol=5)
  for (secdIndx in 1:ncol(pACAT.add_secondOrder.M)) {
    secondG_k2 <- as.matrix(temp.G_k2.M)[ , ((secdIndx-1)*ncol(G)+1):(secdIndx*ncol(G))]
    pACAT.add_secondOrder.M[ , secdIndx] <- Get.p.modifiedGLM(G=secondG_k2, obj_nullmodel_glm=obj_nullmodel_identity)
  }
  pACAT.add_HMM <- matrix(Get.p.modifiedGLM(G=G_k3, obj_nullmodel_glm=obj_nullmodel_identity), ncol=1)
 
  # single variant Cauchy test summary statistics
  W.single.MK1 <- -log10(p.add)-apply(-log10(p.add_k),1,max)
  W.single.MK2 <- (-log10(p.add)-apply(-log10(p.add_k),1,median))*(W.single.MK1>=0)
  W.single.MK3 <- -log10(p.add)-(-log10(p.add_k[,1]))
  W.single.secOrder1 <- -log10(p.add)-apply(-log10(p.add_k2),1,max)
  W.single.secOrder2 <- (-log10(p.add)-apply(-log10(p.add_k2),1,median))*(W.single.secOrder1>=0)
  W.single.secOrder3 <- -log10(p.add)-(-log10(p.add_k2[,1]))
  W.single.HMM1 <- -log10(p.add)-apply(-log10(p.add_k3),1,max)
  W.single.HMM2 <- (-log10(p.add)-apply(-log10(p.add_k3),1,median))*(W.single.HMM1>=0)
  W.single.HMM3 <- -log10(p.add)-(-log10(p.add_k3[,1]))
  W.single.STAAR1 <- -log10(pSTAAR.add)-apply(-log10(pSTAAR.add_k),1,max)
  W.single.STAAR2 <- (-log10(pSTAAR.add)-apply(-log10(pSTAAR.add_k),1,median))*(W.single.STAAR1>=0)
  # W.single.STAAR3 <- -log10(pSTAAR.add_identity)-apply(-log10(pSTAAR.add_k_identity),1,max)
  # W.single.STAAR4 <- (-log10(pSTAAR.add_identity)-apply(-log10(pSTAAR.add_k_identity),1,median))*(W.single.STAAR3>=0)
  W.multiple.secOrder1 <- -log10(p.add)-apply(-log10(pACAT.add_secondOrder.M),1,max)
  W.multiple.secOrder2 <- (-log10(p.add)-apply(-log10(pACAT.add_secondOrder.M),1,median))*(W.multiple.secOrder1>=0)
  W.multiple.secOrder3 <- -log10(pSTAAR.add)-apply(-log10(pSecondOrderSTAAR.add_k),1,max)
  W.multiple.secOrder4 <- (-log10(pSTAAR.add)-apply(-log10(pSecondOrderSTAAR.add_k),1,median))*(W.multiple.secOrder3>=0)

  W.single <- cbind(W.single.MK1, W.single.MK2, W.single.MK3, W.single.secOrder1, W.single.secOrder2, W.single.secOrder3,
                    W.single.HMM1, W.single.HMM2, W.single.HMM3, W.single.STAAR1, W.single.STAAR2,
                    W.multiple.secOrder1, W.multiple.secOrder2, W.multiple.secOrder3, W.multiple.secOrder4); W.single[is.na(W.single)] <- 0
  # knockoff-based inference for Z scores with/out sample relatedness adjustment
  W.ZScores <- cbind(W_summaryStat_score, W_summaryStat_STAAR, W_summaryStat_STAARAdjust,
                     W_ghostScore, W_ghostScore.STAAR); W.ZScores[is.na(W.ZScores)]<-0
  
  # filter single variants
  singleResult <- c()
  for (fdr in seq(0,0.2,by=0.01)) {
    thres.single <- c(MK.threshold(-log10(p.add),-log10(p.add_k),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.MK3, fdr=fdr,offset=1),
                      MK.threshold(-log10(p.add),-log10(p.add_k2),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k2),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.secOrder3, fdr=fdr,offset=1),
                      MK.threshold(-log10(p.add),-log10(p.add_k3),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k3),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.HMM3, fdr=fdr,offset=1),
                      MK.threshold(-log10(pSTAAR.add),-log10(pSTAAR.add_k),fdr=fdr,method='max'),
                      MK.threshold(-log10(pSTAAR.add),-log10(pSTAAR.add_k),fdr=fdr,method='median'),
                      MK.threshold(-log10(p.add),-log10(pACAT.add_secondOrder.M),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(pACAT.add_secondOrder.M),fdr=fdr,method='median'),
                      MK.threshold(-log10(pSTAAR.add),-log10(pSecondOrderSTAAR.add_k),fdr=fdr,method='max'),
                      MK.threshold(-log10(pSTAAR.add),-log10(pSecondOrderSTAAR.add_k),fdr=fdr,method='median'))
    # MK.threshold(-log10(pSTAAR.add_identity),-log10(pSTAAR.add_k_identity),fdr=fdr,method='max'),
    # MK.threshold(-log10(pSTAAR.add_identity),-log10(pSTAAR.add_k_identity),fdr=fdr,method='median'))
    Get.singleSelect <- function(x) {W.single[ , x] >= thres.single[x]}
    singleSelected <- sapply(1:ncol(W.single), Get.singleSelect)
    singleSelected[is.na(singleSelected)] <- F
    
    Get.ZScoreSelect <- function(x){W.ZScores[,x] <= fdr}
    ZScoreSelected <- sapply(1:ncol(W.ZScores), Get.ZScoreSelect)
    ZScoreSelected[is.na(ZScoreSelected)] <- F
    
    singleSignal <- (beta!=0)
    
    singlePower <- function(x) {sum(singleSignal&x)/max(1,sum(singleSignal))}
    singleFdp <- function(x) {sum((!singleSignal)&x)/max(1,sum(x))}
    
    singleResult <- rbind(singleResult, c(apply(cbind(singleSelected,ZScoreSelected),2,singlePower), 
                                          apply(cbind(singleSelected,ZScoreSelected),2,singleFdp)))
    #singleResult <- rbind(singleResult, c(apply(singleSelected,2,singlePower), apply(singleSelected,2,singleFdp)))
  }
  if (replicateIndx==1) {singleResult.fdr <- singleResult} else {singleResult.fdr <- singleResult.fdr + singleResult}
  print(replicateIndx)
  single.temp.print <- cbind(diseaseSnip, seq(0,0.2,by=0.01), singleResult.fdr/replicateIndx)
  print(single.temp.print)
  
  singleResult.detail<-rbind(singleResult.detail, cbind(replicateIndx,h2,N.causal,seq(0,0.2,by=0.01),singleResult))
  colnames(singleResult.detail)<-c('replicate','h2','N.causal','fdr',
                                   'MSK_max_single_power','MSK_median_single_power','SK_single_power',
                                   'NormalK_max_single_power','NormalK_median_single_power', 'NormalSK_single_power',
                                   'HMMK_max_single_power','HMMK_median_single_power', 'HMMSK_single_power',
                                   'STAAR_max_single_power', 'STAAR_median_single_power',
                                   'NormalK_max_multipleGLM_power','NormalK_median_multipleGLM_power', 
                                   'NormalK_max_multipleGLMM_power','NormalK_median_multipleGLMM_power', 
                                   'SummaryStat_score_power', 'SummaryStat_STAAR_power', 'ProposedAdjust_STAAR_power',
                                   'Ghostknockoff_score_power', 'Ghostknockoff_STAAR_power', 
                                   # 'STAAR_identity_max_single_power', 'STAAR_identity_median_single_power',
                                   'MSK_max_single_fdp','MSK_median_single_fdp','SK_single_fdp',
                                   'NormalK_max_single_fdp','NormalK_median_single_fdp', 'NormalSK_single_fdp',
                                   'HMMK_max_single_fdp','HMMK_median_single_fdp', 'HMMSK_single_fdp',
                                   'STAAR_max_single_fdp', 'STAAR_median_single_fdp',
                                   'NormalK_max_multipleGLM_fdp','NormalK_median_multipleGLM_fdp', 
                                   'NormalK_max_multipleGLMM_fdp','NormalK_median_multipleGLMM_fdp', 
                                   'SummaryStat_score_fdp', 'SummaryStat_STAAR_fdp', 'ProposedAdjust_STAAR_fdp',
                                   'Ghostknockoff_score_fdp', 'Ghostknockoff_STAAR_fdp')
  singleResult.detail <- round(singleResult.detail,digits=3)
  write.table(singleResult.detail,paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/',location,'/singleResult_knockoff_',out_type, '_',replicateIndx,'_',diseaseSnip,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}
