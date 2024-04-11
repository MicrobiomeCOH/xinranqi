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
library(aod)
library(lmtest)

#setwd('/scratch/users/xinranqi/Project1/CCF')
#source('KnockoffScreen_Nov5.R')
#source('ESGWAS.R')
#setwd('/Users/xqi/Documents/Week2_01112022/family_simulation')
#source('/Users/xqi/Documents/Week6_02132022/KnockoffScreen_Nov5.R')
source('/oak/stanford/groups/zihuai/XinranQi/Project1/simulate10PedigreeRecombine0.R')
source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
library(KnockoffScreen)

data("SKAT.haplotypes")
SKATHaplotypes <- SKAT.haplotypes$Haplotype[ , SKAT.haplotypes$SNPInfo$FREQ1>0.001]
positionChromosome <- SKAT.haplotypes$SNPInfo$CHROM_POS[SKAT.haplotypes$SNPInfo$FREQ1>0.001]

set.seed(123)

result.all <- result.detail <- singleResult.detail <- combinedResult.detail <- c()

for (replicateIndx in 1:1000){
  # select some snips to be included in the simulation
  selectedSNPsIndx <- sort(sample(ncol(SKATHaplotypes), size=numOfSelectedSnips, replace=FALSE))
  chromosomes <- SKATHaplotypes[sample(1:nrow(SKATHaplotypes), size=2*(numIndependentDiseased+numIndependentControl), replace=TRUE), selectedSNPsIndx]
  # convert to genotypes (coding 0,1,2) #generate genotype from haplotype
  G <- chromosomes[1:(numIndependentDiseased+numIndependentControl), ] + chromosomes[(numIndependentDiseased+numIndependentControl+1):(2*(numIndependentDiseased+numIndependentControl)), ] 

  MAF <- apply(G,2,mean)/2
  G[,MAF>0.5] <- 2-G[,MAF>0.5]
  MAF <- apply(G,2,mean)/2
  MAC <- apply(G,2,sum)
  sd.G <- apply(G,2,sd)
  #filtering based on MAC
  G.index <- MAF>=0.01 & sd.G>0
  G <- G[,G.index]; MAF <- MAF[G.index]; MAC <- MAC[G.index]
  #pos <- positionChromosome[selectedSNPsIndx]
  # case control family study
  #MAF <- apply(G,2,mean)/2
  # clustering so that causal variants are in separate SNP blocks
  cor.X <- cor(G)
  #cor.X <- sparse.cor(Matrix(G[,MAF>0.05]))$cor 
  Sigma.distance <- as.dist(1-abs(cor.X))
  fit <- hclust(Sigma.distance, method="single")
  corr_max <- 0.75
  clusters <- cutree(fit, h=1-corr_max)
  # sample SNPs blocks and causal variants
  blockIndx <- sample(unique(clusters), size=N.causal, replace=FALSE)
  SNPsIndx <- sort(sapply(blockIndx, function(x) {ifelse(sum(clusters==x)==1, which(clusters==x), sample(which(clusters==x), size=1))}))
  
  # priority 1: causal variants
  cluster.indexPriority <- SNPsIndx[sapply(unique(clusters[SNPsIndx]), function(x) {min(which(clusters[SNPsIndx]==x))})]
  # cluster.indexSecond <- match(unique(clusters), clusters)
  cluster.indexSecond <- match(setdiff(unique(clusters), unique(clusters[SNPsIndx])), clusters)
  cluster.index <- sort(c(cluster.indexPriority, cluster.indexSecond))
  
  # causal variant index before clustering
  causal.index <- rep(FALSE, times=ncol(G)); causal.index[SNPsIndx] <- TRUE
  causal.index[sample(which(causal.index>0), N.causal/2, replace=FALSE)] <- -1
  
  # keep representative SNPs from SNPs blocks
  G <- G[ , cluster.index, drop=F]
  pos <- ((positionChromosome[selectedSNPsIndx])[G.index])[cluster.index]
  # causal variant index after clustering
  causal.index <- causal.index[cluster.index]
  # regression coefficients
  beta <- sqrt(h2/N.causal/apply(G,2,var))*causal.index
  
  # baseline/common disease rate
  prevalence <- 0.1
  b0_pheno <- log(prevalence/(1-prevalence)) # for binary trait
  
  # confounding covariate
  X1 <- rnorm(nrow(G),0,1)
  mu <- X1+G%*%beta+rnorm(length(X1),0,1)
  mu <- mu-mean(mu)
  #y <- sapply(mu, function(x) {rbinom(1,1,(exp(x+b0_pheno)/(1+exp(x+b0_pheno))))})
  y <- as.matrix(rbinom((numIndependentDiseased+numIndependentControl),1, prob=as.numeric(exp(mu+b0_pheno)/(1+exp(mu+b0_pheno)))))

  # random sample from it for equal number of cases/controls
  caseIndx <- sample(which(y==1), 5000, replace=FALSE)
  controlIndx <- sample(which(y==0), 5000, replace=FALSE)
  G <- G[c(caseIndx, controlIndx), ]
  X1 <- X1[c(caseIndx, controlIndx)]
  y <- y[c(caseIndx, controlIndx)]
  
  MAF <- apply(G[y==0,],2,mean)/2
  G[,MAF>0.5] <- 2-G[,MAF>0.5]
  MAF <- apply(G[y==0,],2,mean)/2

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
  #X_k <- apply(G_k,1,function(x) x%*%window.matrix)
  X_k2 <- as.matrix(G_k2%*%window.matrix)
  X_k3 <- as.matrix(G_k3%*%window.matrix)
  
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
  score_k<-as.data.frame(lapply(G_k,function(x) t(as.matrix(x))%*%Y.res))
  #score_k<-apply(G_k,1,function(x) t(x)%*%Y.res);re.score_k<-lapply(1:dim(G_k)[1],function(s) t(t(G_k[s,,])%*%re.Y.res))
  re.score_k<-lapply(G_k,function(s) t(t(as.matrix(s))%*%re.Y.res))

  #all Cauchy test
  MAF<-apply(G,2,mean)/2
  MAC<-apply(G,2,sum)
  weight.beta<-dbeta(apply(G,2,mean)/2,1,25)
  weight.matrix<-cbind(MAC<10,(MAF<0.01&MAC>10)*weight.beta,(MAF>=0.01)*weight.beta)
  colnames(weight.matrix)<-c('MAC<5','MAF<0.01&MAC>5&Beta','MAF>=0.01Beta')
  # the updated KS.test remove Get.p() steps
  p.add<-as.matrix(Get.p(G,result.prelim))
  p.add_k<-as.matrix(as.data.frame(lapply(G_k,function(x) Get.p(as.matrix(x),result.prelim=result.prelim))))
  #p.single<-KS.fit$p.single
  #p.single_k<-KS.fit$p.single_k
  #p.add<-as.matrix(p.single); p.add_k<-p.single_k
  KS.fit<-KS.test(temp.G=G,temp.G_k=G_k,p.rare=p.add,p.rare_k=p.add_k,result.prelim=result.prelim,window.matrix=window.matrix0,weight.matrix=weight.matrix)
  
  # SummaryStat knockoff - score test
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
  
  p.A <- KS.fit$p.KS
  p.A_k <- KS.fit$p.KS_k
  p.A.combined <- rbind(p.A, p.add)
  p.A_k.combined <- rbind(p.A_k, p.add_k)
  
  #burden
  p.burden <- Get.p.base(X,result.prelim)
  p.burden_k <- sapply(X_k,Get.p.base,result.prelim=result.prelim)
  p.burden.combined <- rbind(p.burden, p.add)
  p.burden_k.combined <- rbind(p.burden_k, p.add_k)
  
  # dispersion
  p.dispersion <- Get.p.moment(as.vector(t(score^2)%*%window.matrix^2),re.score^2%*%window.matrix^2)
  p.dispersion_k <- sapply(1:ncol(score_k),function(s){Get.p.moment(as.vector(t(score_k[,s]^2)%*%window.matrix^2),re.score_k[[s]]^2%*%window.matrix^2)})
  p.dispersion.combined <- rbind(p.dispersion, p.add)
  p.dispersion_k.combined <- rbind(p.dispersion_k, p.add_k)

  # Knockoff statistics
  # window-based summary statistics only
  W1<--log10(p.burden)-apply(-log10(p.burden_k),1,max)
  W2<--log10(p.dispersion)-apply(-log10(p.dispersion_k),1,max)
  W3<--log10(p.A)-apply(-log10(p.A_k),1,max)
  W4<-(-log10(p.burden)-apply(-log10(p.burden_k),1,median))*(W1>=0)
  W5<-(-log10(p.dispersion)-apply(-log10(p.dispersion_k),1,median))*(W2>=0)
  W6<-(-log10(p.A)-apply(-log10(p.A_k),1,median))*(W3>=0)
  W7<--log10(p.burden)-(-log10(p.burden_k[,1]))
  W8<--log10(p.dispersion)-(-log10(p.dispersion_k[,1]))
  W9<--log10(p.A)-(-log10(p.A_k[,1]))
  
  # combined summary statistics
  W.combine1 <- -log10(p.burden.combined)-apply(-log10(p.burden_k.combined),1,max)
  W.combine2 <- -log10(p.dispersion.combined)-apply(-log10(p.dispersion_k.combined),1,max)
  W.combine3 <- -log10(p.A.combined)-apply(-log10(p.A_k.combined),1,max)
  W.combine4 <- (-log10(p.burden.combined)-apply(-log10(p.burden_k.combined),1,median))*(W.combine1>=0)
  W.combine5 <- (-log10(p.dispersion.combined)-apply(-log10(p.dispersion_k.combined),1,median))*(W.combine2>=0)
  W.combine6 <- (-log10(p.A.combined)-apply(-log10(p.A_k.combined),1,median))*(W.combine3>=0)
  W.combine7 <- -log10(p.burden.combined)-(-log10(p.burden_k.combined[,1]))
  W.combine8 <- -log10(p.dispersion.combined)-(-log10(p.dispersion_k.combined[,1]))
  W.combine9 <- -log10(p.A.combined)-(-log10(p.A_k.combined[,1]))
  
  p.add_k2 <- as.matrix(Get.p(G_k2,result.prelim=result.prelim))
  p.add_k3 <- as.matrix(Get.p(G_k3,result.prelim=result.prelim))
  
  p.burden_k2 <- Get.p.base(X_k2,result.prelim)
  p.burden_k2.combined <- rbind(p.burden_k2, p.add_k2)
  W10 <- -log10(p.burden)-(-log10(p.burden_k2))
  W.combine10 <- -log10(p.burden.combined)-(-log10(p.burden_k2.combined))
  
  score_k2<-t(G_k2)%*%Y.res;re.score_k2<-t(t(G_k2)%*%re.Y.res)
  p.dispersion_k2 <- Get.p.moment(as.vector(t(score_k2^2)%*%window.matrix^2),re.score_k2^2%*%window.matrix^2)
  p.dispersion_k2.combined <- rbind(p.dispersion_k2,p.add_k2)
  W11 <- -log10(p.dispersion)-(-log10(p.dispersion_k2))
  W.combine11 <- -log10(p.dispersion.combined)-(-log10(p.dispersion_k2.combined))
  
  p.burden_k3 <- Get.p.base(X_k3,result.prelim)
  p.burden_k3.combined <- rbind(p.burden_k3,p.add_k3)
  W12 <- -log10(p.burden)-(-log10(p.burden_k3))
  W.combine12 <- -log10(p.burden.combined)-(-log10(p.burden_k3.combined))
  
  score_k3<-t(G_k3)%*%Y.res;re.score_k3<-t(t(G_k3)%*%re.Y.res)
  p.dispersion_k3 <- Get.p.moment(as.vector(t(score_k3^2)%*%window.matrix^2),re.score_k3^2%*%window.matrix^2)
  p.dispersion_k3.combined <- rbind(p.dispersion_k3,p.add_k3)
  W13 <- -log10(p.dispersion)-(-log10(p.dispersion_k3))
  W.combine13 <- -log10(p.dispersion.combined)-(-log10(p.dispersion_k3.combined))

  phenoSTAAR <- data.frame(Y=y,X1=X1,id=1:length(y))
  # ACAT without adjustment for kinship
  # window based using ACAT
  obj <- NULL_Model(as.vector(y),Z=as.matrix(X1))
  pValuesACAT <- sapply(1:ncol(window.matrix), function(x) {ACAT_V(G[ , window.matrix[ , x]!=0],obj)})
  # SCIT knockoff, ACAT & SCIP knockoff, ACAT
  # apply STAAR to knockoff counterparts with identity kinship matrix GLM model & follow the same knockoff filtering procedure 
  pValuesACAT_k <- matrix(NA, nrow=length(pValuesACAT), ncol=length(G_k))
  for (repIndx in 1:length(G_k)) {
  #for (repIndx in 1:(dim(G_k)[1])) {
    pValuesACAT_k[ , repIndx] <- sapply(1:length(pValuesACAT), function(x) {ACAT_V(as.matrix(G_k[[repIndx]])[ , window.matrix[,x]!=0],obj)})
    #pValuesACAT_k[ , repIndx] <- sapply(1:length(pValuesACAT), function(x) {ACAT_V(G_k[repIndx, , window.matrix[ , x]!=0],obj)})
  }
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
  pValuesACAT.combined <- rbind(matrix(pValuesACAT, ncol=1), pACAT.add)
  pValuesACAT_k.combined <- rbind(pValuesACAT_k, pACAT.add_k)
  
  # correct for 0 p-value, otherwise Zscore_0 is Inf and GK.stat T_0/T_k become -Inf
  #pACAT.add[pACAT.add==0] <- min(c(pACAT.add[pACAT.add>0], 1e-100))
  
  # SummaryStat knockoff - mixed model score test
  n<-length(obj_nullmodel_identity$residuals)
  ZscoreAdd_0<-qnorm(pACAT.add/2,lower.tail=F)*sign(cov(as.matrix(G),obj_nullmodel_identity$residuals)) #convert p-values to z-scores
  ESSTAAR.stat<-ES.fit(ZscoreAdd_0,n,gamma=1,fit.prelim=fit.prelim.M) #fit summary stat based knockoff
  MKSTAAR.stat<-MK.statistic(ESSTAAR.stat$T_0,ESSTAAR.stat$T_k) #compute knockoff statistics
  W_summaryStat_STAAR<-MK.q.byStat(MKSTAAR.stat[,'kappa'],MKSTAAR.stat[,'tau'],M=5,Rej.Bound=10000) #calculate q value

  # Ghostknockoff mixed model score test
  GK.stat.ACAT <- GhostKnockoff.fit(Zscore_0=ZscoreAdd_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  #GK.stat.ACAT <- GhostKnockoff.fit(ZscoreAdd_0,n,fit.prelim=fit.prelim.Ghost,gamma=1,weight.study=NULL)
  GK.filter.ACAT <- MK.statistic(GK.stat.ACAT$T_0[[1]],GK.stat.ACAT$T_k[[1]])
  #GK.filter.ACAT <- GhostKnockoff.filter(GK.stat.ACAT$T_0,GK.stat.ACAT$T_k)
  #W_ghostScore.ACAT <- rep(1,nrow(GK.stat.ACAT$T_0))
  W_ghostScore.ACAT <- MK.q.byStat(GK.filter.ACAT[,'kappa'],GK.filter.ACAT[,'tau'],5,Rej.Bound=10000)
  #W_ghostScore.ACAT <- GK.filter.ACAT$q
  
  # second order knockoff, ACAT
  pValuesACAT_secondOrder <- sapply(1:length(pValuesACAT), 
                                    function(x) {ifelse(sum(window.matrix[ , x]!=0) == 1, 
                                                        Get.p.modifiedGLM(G=G_k2[ , window.matrix[ , x]!=0], obj_nullmodel_glm=obj_nullmodel_identity),
                                                        ACAT_V(G_k2[ , window.matrix[ , x]!=0],obj))})
  # multiple second order, ACAT
  pValuesACAT_secondOrder.M <- matrix(NA, nrow=length(pValuesACAT_secondOrder), ncol=5)
  for (secdIndx in 1:ncol(pValuesACAT_secondOrder.M)) {
    secondG_k2 <- as.matrix(temp.G_k2.M)[ , ((secdIndx-1)*ncol(G)+1):(secdIndx*ncol(G))]
    pValuesACAT_secondOrder.M[ , secdIndx] <- sapply(1:length(pValuesACAT), 
                                                     function(x) {ifelse(sum(window.matrix[ , x]!=0) == 1, 
                                                                         Get.p.modifiedGLM(G=secondG_k2[ , window.matrix[ , x]!=0], obj_nullmodel_glm=obj_nullmodel_identity),
                                                                         ACAT_V(secondG_k2[ , window.matrix[ , x]!=0],obj))})
  }
  # HMM knockoff, ACAT
  pValuesACAT_HMM <- sapply(1:length(pValuesACAT), 
                            function(x) {ifelse(sum(window.matrix[ , x]!=0) == 1, 
                                                Get.p.modifiedGLM(G=G_k3[ , window.matrix[ , x]!=0], obj_nullmodel_glm=obj_nullmodel_identity),
                                                ACAT_V(G_k3[ , window.matrix[ , x]!=0],obj))})
  pACAT.add_secondOrder <- matrix(Get.p.modifiedGLM(G=G_k2, obj_nullmodel_glm=obj_nullmodel_identity), ncol=1)
  # multiple second order, GLM single variant
  pACAT.add_secondOrder.M <- matrix(NA, nrow=length(pACAT.add_secondOrder), ncol=5)
  for (secdIndx in 1:ncol(pACAT.add_secondOrder.M)) {
    secondG_k2 <- as.matrix(temp.G_k2.M)[ , ((secdIndx-1)*ncol(G)+1):(secdIndx*ncol(G))]
    pACAT.add_secondOrder.M[ , secdIndx] <- Get.p.modifiedGLM(G=secondG_k2, obj_nullmodel_glm=obj_nullmodel_identity)
  }
  pACAT.add_HMM <- matrix(Get.p.modifiedGLM(G=G_k3, obj_nullmodel_glm=obj_nullmodel_identity), ncol=1)
  pValuesACAT_k_secondOrderComb <- rbind(matrix(pValuesACAT_secondOrder, ncol=1), pACAT.add_secondOrder)
  pValuesACAT_k_secondOrderComb.M <- rbind(pValuesACAT_secondOrder.M, pACAT.add_secondOrder.M)
  pValuesACAT_k_HMMComb <- rbind(matrix(pValuesACAT_HMM, ncol=1), pACAT.add_HMM)
  
  # GLM + identity kinship matrix
  W16 <- -log10(pValuesACAT)-apply(-log10(pValuesACAT_k),1,max)
  # SCIT knockoffs, ACAT
  W17 <- (-log10(pValuesACAT)-apply(-log10(pValuesACAT_k),1,median))*(W16>=0)
  # SCIP knockoff, ACAT
  W18 <- -log10(pValuesACAT)-(-log10(pValuesACAT_k)[ , 1])
  W.combine16 <- -log10(pValuesACAT.combined)-apply(-log10(pValuesACAT_k.combined),1,max)
  W.combine17 <- (-log10(pValuesACAT.combined)-apply(-log10(pValuesACAT_k.combined),1,median))*(W.combine16>=0)
  W.combine18 <- -log10(pValuesACAT.combined)-(-log10(pValuesACAT_k.combined)[ , 1])
  
  # window-based & combined ACAT score test
  # second order knockoff, ACAT
  W19 <- -log10(pValuesACAT)-(-log10(pValuesACAT_secondOrder))
  # HMM knockoff, ACAT
  W20 <- -log10(pValuesACAT)-(-log10(pValuesACAT_HMM))
  W.combine19 <- -log10(pValuesACAT.combined)-(-log10(pValuesACAT_k_secondOrderComb))
  W.combine20 <- -log10(pValuesACAT.combined)-(-log10(pValuesACAT_k_HMMComb))
  # multiple second order, ACAT
  W21 <- -log10(pValuesACAT)-apply(-log10(pValuesACAT_secondOrder.M),1,max)
  W22 <- (-log10(pValuesACAT)-apply(-log10(pValuesACAT_secondOrder.M),1,median))*(W21>=0)
  W.combine21 <- -log10(pValuesACAT.combined)-apply(-log10(pValuesACAT_k_secondOrderComb.M),1,max)
  W.combine22 <- (-log10(pValuesACAT.combined)-apply(-log10(pValuesACAT_k_secondOrderComb.M),1,median))*(W.combine21>=0)
  
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
  #W.single.STAAR1 <- -log10(pSTAAR.add)-apply(-log10(pSTAAR.add_k),1,max)
  #W.single.STAAR2 <- (-log10(pSTAAR.add)-apply(-log10(pSTAAR.add_k),1,median))*(W.single.STAAR1>=0)
  W.single.STAAR3 <- -log10(pACAT.add)-apply(-log10(pACAT.add_k),1,max)
  W.single.STAAR4 <- (-log10(pACAT.add)-apply(-log10(pACAT.add_k),1,median))*(W.single.STAAR3>=0)
  W.multiple.secOrder1 <- -log10(p.add)-apply(-log10(pACAT.add_secondOrder.M),1,max)
  W.multiple.secOrder2 <- (-log10(p.add)-apply(-log10(pACAT.add_secondOrder.M),1,median))*(W.multiple.secOrder1>=0)

  W<-cbind(W1,W2,W3,W4,W5,W6,W7,W8,W9,W10,W11,W12,W13,W16,W17,W18,W19,W20,W21,W22); W[is.na(W)]<-0
  W.single <- cbind(W.single.MK1, W.single.MK2, W.single.MK3, W.single.secOrder1, W.single.secOrder2, W.single.secOrder3,
                    W.single.HMM1, W.single.HMM2, W.single.HMM3, W.single.STAAR3, W.single.STAAR4,
                    W.multiple.secOrder1, W.multiple.secOrder2); W.single[is.na(W.single)] <- 0
  
  
  # add LRT & Wald test p-values for ghost knockoffs
  logitNull <- glm(as.factor(Y)~X1, data=phenoSTAAR, family="binomial")
  waldP.add <- LRTp.add <- numeric(ncol(G))
  for (genoIndx in 1:ncol(G)) {
    phenoSTAARGeno <- cbind.data.frame(phenoSTAAR, as.factor(G[,genoIndx]))
    colnames(phenoSTAARGeno) <- c(colnames(phenoSTAAR), "geno")
    # test for SNPs
    logitFull <- glm(as.factor(Y)~X1+geno, data=phenoSTAARGeno, family="binomial")
    # Wald test
    waldTemp <- ifelse(length(coef(logitFull))==3, summary(logitFull)$coefficients[3,"Pr(>|z|)"], 
                       wald.test(b=coef(logitFull), Sigma=vcov(logitFull), Terms=3:length(coef(logitFull)))$result$chi2["P"])
    waldP.add[genoIndx] <- if (waldTemp == 0) {
      Lvec <- if (length(coef(logitFull)) == 4) {
        matrix(c(0,0,1,0,0,0,0,1), nrow=2, ncol=4, byrow=TRUE)
      } else {
        matrix(c(0,0,1), nrow=1)
      }
      betaVec <- as.matrix(coef(logitFull))
      exp(pchisq(t(Lvec%*%betaVec) %*% solve(Lvec%*%vcov(logitFull)%*%t(Lvec)) %*% (Lvec%*%betaVec), 
                 df=2, lower.tail=FALSE, log.p=TRUE))
    } else {
      waldTemp
    }
    LRTp.add[genoIndx] <- tail(lrtest(logitFull, logitNull)$`Pr(>Chisq)`, n=1)
  }
  
  # Ghostknockoff from Wald test
  n <- length(phenoSTAAR$Y)
  WaldY.res <- logitNull$data$Y - logitNull$fitted.values
  WaldZscore_0 <- qnorm(waldP.add/2,lower.tail=F)*sign(cov(as.matrix(G),WaldY.res)) #convert p-values to z-scores
  # 0 p-value convert to Inf z score, replace with biggest absolute value
  #WaldZscore_0 <- as.matrix(sapply(WaldZscore_0, function(x) {ifelse(is.infinite(x), max(abs(WaldZscore_0[!is.infinite(WaldZscore_0)]))*sign(x), x)}))
  WaldGK.stat <- GhostKnockoff.fit(Zscore_0=WaldZscore_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  WaldGK.filter <- MK.statistic(WaldGK.stat$T_0[[1]],WaldGK.stat$T_k[[1]])
  #WaldW_ghostScore <- rep(1,nrow(WaldGK.stat$T_0))
  WaldW_ghostScore <- MK.q.byStat(WaldGK.filter[,'kappa'],WaldGK.filter[,'tau'],5,Rej.Bound=10000)
  
  # Ghostknockoff from likelihood ratio test
  LRTZscore_0 <- qnorm(LRTp.add/2,lower.tail=F)*sign(cov(as.matrix(G),WaldY.res)) #convert p-values to z-scores
  # 0 p-value convert to Inf z score, replace with biggest absolute value
  LRTZscore_0 <- as.matrix(sapply(LRTZscore_0, function(x) {ifelse(is.infinite(x), max(abs(LRTZscore_0[!is.infinite(LRTZscore_0)]))*sign(x), x)}))
  LRTGK.stat <- GhostKnockoff.fit(Zscore_0=LRTZscore_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  LRTGK.filter <- MK.statistic(LRTGK.stat$T_0[[1]],LRTGK.stat$T_k[[1]])
  #LRTW_ghostScore <- rep(1,nrow(LRTGK.stat$T_0))
  LRTW_ghostScore <- MK.q.byStat(LRTGK.filter[,'kappa'],LRTGK.filter[,'tau'],5,Rej.Bound=10000)
  
  
  # knockoff-based inference for Z scores with/out sample relatedness adjustment
  W.ZScores <- cbind(W_summaryStat_score, W_summaryStat_STAAR, W_ghostScore, W_ghostScore.ACAT,
                     WaldW_ghostScore, LRTW_ghostScore); W.ZScores[is.na(W.ZScores)]<-0
  
  W.combined <- cbind(W.combine1, W.combine2, W.combine3, W.combine4, W.combine5, W.combine6, W.combine7, W.combine8, W.combine9,
                      W.combine10, W.combine11, W.combine12, W.combine13, 
                      W.combine16, W.combine17, W.combine18, W.combine19, W.combine20,
                      W.combine21, W.combine22); W.combined[is.na(W.combined)] <- 0
  
  # single variable knockoff-based inference methods
  singleResult <- c()
  # add model-X knockoffs to determine the reason
  modelXMatrix <- cbind.data.frame(X1, G)
  modelXCopy <- cbind(rnorm(length(X1)), G_k2)
  modelXW <- stat.glmnet_coefdiff(X=modelXMatrix, X_k=modelXCopy, y=y, nfolds=10, family="binomial", penalty.factor=rep(c(0, rep(1, times=ncol(G))), times=2))

  for (fdr in seq(0,0.2,by=0.01)) {
    #modelXresult <- knockoff.filter(modelXMatrix, y, fdr=fdr)
    #modelXselect <- if (length(modelXresult$selected)>0) {as.numeric(names(modelXresult$selected)[names(modelXresult$selected)!="X1"])} else {NULL}
    #modelXSelected <- rep(FALSE, times=length(beta)); modelXSelected[modelXselect] <- TRUE
    modelXThres <- knockoff.threshold(modelXW, fdr=fdr, offset=1)
    modelXSelected <- (modelXW>=modelXThres)[2:length(modelXW)]
    
    thres.single <- c(MK.threshold(-log10(p.add),-log10(p.add_k),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.MK3, fdr=fdr,offset=1),
                      MK.threshold(-log10(p.add),-log10(p.add_k2),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k2),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.secOrder3, fdr=fdr,offset=1),
                      MK.threshold(-log10(p.add),-log10(p.add_k3),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(p.add_k3),fdr=fdr,method='median'),
                      knockoff.threshold(W.single.HMM3, fdr=fdr,offset=1),
                      #MK.threshold(-log10(pSTAAR.add),-log10(pSTAAR.add_k),fdr=fdr,method='max'),
                      #MK.threshold(-log10(pSTAAR.add),-log10(pSTAAR.add_k),fdr=fdr,method='median'))
                      MK.threshold(-log10(pACAT.add),-log10(pACAT.add_k),fdr=fdr,method='max'),
                      MK.threshold(-log10(pACAT.add),-log10(pACAT.add_k),fdr=fdr,method='median'),
                      MK.threshold(-log10(p.add),-log10(pACAT.add_secondOrder.M),fdr=fdr,method='max'),
                      MK.threshold(-log10(p.add),-log10(pACAT.add_secondOrder.M),fdr=fdr,method='median'))
    Get.singleSelect <- function(x) {W.single[ , x] >= thres.single[x]}
    singleSelected <- sapply(1:ncol(W.single), Get.singleSelect)
    singleSelected[is.na(singleSelected)] <- F
    
    Get.ZScoreSelect <- function(x){W.ZScores[,x] <= fdr}
    ZScoreSelected <- sapply(1:ncol(W.ZScores), Get.ZScoreSelect)
    ZScoreSelected[is.na(ZScoreSelected)] <- F
    #min(sapply(1:(sum(beta>0)-1), function(x){pos[which(beta>0)[x+1]]-pos[which(beta>0)[x]]}))
    #bmatrix(round(cor(G[,which(causal.index)]), digits=3))
    
    singleSignal <- (beta!=0)
    
    singlePower <- function(x) {sum(singleSignal&x)/max(1,sum(singleSignal))}
    singleFdp <- function(x) {sum((!singleSignal)&x)/max(1,sum(x))}
    
    singleResult <- rbind(singleResult, c(apply(cbind(singleSelected,ZScoreSelected,modelXSelected),2,singlePower), 
                                          apply(cbind(singleSelected,ZScoreSelected,modelXSelected),2,singleFdp)))
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
                                   #'STAAR_max_single_power', 'STAAR_median_single_power',
                                   'STAAR_identity_max_single_power', 'STAAR_identity_median_single_power',
                                   'NormalK_max_multipleGLM_power','NormalK_median_multipleGLM_power', 
                                   'SummaryStat_score_power', 'SummaryStat_ACAT_power', 
                                   'Ghostknockoff_score_power', 'Ghostknockoff_ACAT_power', 
                                   'Ghostknockoff_Wald_power', 'Ghostknockoff_LRT_power', 'ModelX_power',
                                   'MSK_max_single_fdp','MSK_median_single_fdp','SK_single_fdp',
                                   'NormalK_max_single_fdp','NormalK_median_single_fdp', 'NormalSK_single_fdp',
                                   'HMMK_max_single_fdp','HMMK_median_single_fdp', 'HMMSK_single_fdp',
                                   #'STAAR_max_single_fdp', 'STAAR_median_single_fdp',
                                   'STAAR_identity_max_single_fdp', 'STAAR_identity_median_single_fdp',
                                   'NormalK_max_multipleGLM_fdp','NormalK_median_multipleGLM_fdp', 
                                   'SummaryStat_score_fdp', 'SummaryStat_ACAT_fdp',
                                   'Ghostknockoff_score_fdp', 'Ghostknockoff_ACAT_fdp',
                                   'Ghostknockoff_Wald_fdp', 'Ghostknockoff_LRT_fdp', 'ModelX_fdp')
  singleResult.detail <- round(singleResult.detail,digits=3)
  write.table(singleResult.detail,paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/', location,'/IndepSingleResult_knockoff_',out_type, '_',replicateIndx,'_',diseaseSnip,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')

  # window-based filtering
  result<-c()
  for (fdr in seq(0,0.2,by=0.01)) {
    thres<-c(MK.threshold(-log10(p.burden),-log10(p.burden_k),fdr=fdr,method='max'),
             MK.threshold(-log10(p.dispersion),-log10(p.dispersion_k),fdr=fdr,method='max'),
             MK.threshold(-log10(p.A),-log10(p.A_k),fdr=fdr,method='max'),
             MK.threshold(-log10(p.burden),-log10(p.burden_k),fdr=fdr,method='median'),
             MK.threshold(-log10(p.dispersion),-log10(p.dispersion_k),fdr=fdr,method='median'),
             MK.threshold(-log10(p.A),-log10(p.A_k),fdr=fdr,method='median'),
             apply(W[,-(c(1:6, 14:20))],2,knockoff.threshold,fdr=fdr,offset=1),
             # use MK.threshold() function to finalize corresponding thresholds
             #MK.threshold(-log10(pValuesSTAAR),-log10(pValuesSTAAR_k),fdr=fdr,method='max'),
             #MK.threshold(-log10(pValuesSTAAR),-log10(pValuesSTAAR_k),fdr=fdr,method='median'),
             MK.threshold(-log10(pValuesACAT),-log10(pValuesACAT_k),fdr=fdr,method='max'),
             MK.threshold(-log10(pValuesACAT),-log10(pValuesACAT_k),fdr=fdr,method='median'),
             apply(W[ , 16:18],2,knockoff.threshold,fdr=fdr,offset=1),
             MK.threshold(-log10(pValuesACAT),-log10(pValuesACAT_secondOrder.M),fdr=fdr,method='max'),
             MK.threshold(-log10(pValuesACAT),-log10(pValuesACAT_secondOrder.M),fdr=fdr,method='median'))
    Get.select<-function(x){W[,x] >= thres[x]}
    selected<-sapply(1:ncol(W),Get.select)
    selected[is.na(selected)]<-F
    signal<-t(beta!=0)%*%window.matrix0!=0
    # signal<-c(t(beta!=0)%*%window.matrix0!=0,(beta!=0)[MAF>0.01])
    
    selected.BB<-as.vector(p.burden<(0.05/ncol(X)))
    selected.BD<-as.vector(p.dispersion<(0.05/ncol(X)))
    selected<-cbind(selected,selected.BB,selected.BD)
    
    #BH FDR control
    selected.BHB<-as.vector(p.adjust(p.burden,method='fdr')<=fdr)
    selected.BHD<-as.vector(p.adjust(p.dispersion,method='fdr')<=fdr)
    selected<-cbind(selected,selected.BHB,selected.BHD)

    power<-function(x){sum(signal&x)/max(1,sum(signal))}
    fdp<-function(x){sum((!signal)&x)/max(1,sum(x))}
    
    result<-rbind(result,c(apply(selected,2,power),apply(selected,2,fdp)))
  }
  if(replicateIndx==1){result.fdr <- result}else{result.fdr <- result.fdr+result}
  print(replicateIndx)
  temp.print<-cbind(diseaseSnip,seq(0,0.2,by=0.01),result.fdr/replicateIndx)
  print(temp.print)
  
  result.detail<-rbind(result.detail,cbind(replicateIndx,h2,N.causal,seq(0,0.2,by=0.01),result))
  colnames(result.detail)<-c('replicate','h2','N.causal','fdr',
                             'MSK_max_Burden_power','MSK_max_Dispersion_power','MSK_max_AllCauchy_power',
                             'MSK_median_Burden_power','MSK_median_Dispersion_power','MSK_median_AllCauchy_power',
                             'SK_Burden_power','SK_Dispersion_power','SK_AllCauchy_power',
                             'NormalK_Burden_power','NormalK_Dispersion_power','HMMK_Burden_power','HMMK_Dispersion_power',
                             #'STAAR_max_power', 'STAAR_median_power',
                             'STAAR_identity_max_power', 'STAAR_identity_median_power', 'STAAR_identity_single__power',
                             'NormalK_ACAT_power', 'HMMK_ACAT_power',
                             'NormalK_multipleACAT_max_power', 'NormalK_multipleACAT_median_power',
                             'Bonferroni_Burden_power','Bonferroni_Dispersion_power',
                             'BH_Burden_power','BH_Dispersion_power', 
                             'MSK_max_Burden_fdp','MSK_max_Dispersion_fdp','MSK_max_AllCauchy_fdp',
                             'MSK_median_Burden_fdp','MSK_median_Dispersion_fdp','MSK_median_AllCauchy_fdp',
                             'SK_Burden_fdp','SK_Dispersion_fdp','SK_AllCauchy_fdp',
                             'NormalK_Burden_fdp','NormalK_Dispersion_fdp','HMMK_Burden_fdp','HMMK_Dispersion_fdp',
                             #'STAAR_max_fdp', 'STAAR_median_fdp',
                             'STAAR_identity_max_fdp', 'STAAR_identity_median_fdp', 'STAAR_identity_single__fdp',
                             'NormalK_ACAT_fdp', 'HMMK_ACAT_fdp',
                             'NormalK_multipleACAT_max_fdp', 'NormalK_multipleACAT_median_fdp',
                             'Bonferroni_Burden_fdp','Bonferroni_Dispersion_fdp',
                             'BH_Burden_fdp','BH_Dispersion_fdp')
  result.detail<-round(result.detail,digits=3)
  write.table(result.detail,paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/', location,'/indepResult_knockoff_',out_type, '_',replicateIndx,'_',diseaseSnip,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')
  
  # window-based + added p-values filtering
  combinedResult<-c()
  for (fdr in seq(0,0.2,by=0.01)) {
    thres.combine <- c(MK.threshold(-log10(p.burden.combined),-log10(p.burden_k.combined),fdr=fdr,method='max'),
                       MK.threshold(-log10(p.dispersion.combined),-log10(p.dispersion_k.combined),fdr=fdr,method='max'),
                       MK.threshold(-log10(p.A.combined),-log10(p.A_k.combined),fdr=fdr,method='max'),
                       MK.threshold(-log10(p.burden.combined),-log10(p.burden_k.combined),fdr=fdr,method='median'),
                       MK.threshold(-log10(p.dispersion.combined),-log10(p.dispersion_k.combined),fdr=fdr,method='median'),
                       MK.threshold(-log10(p.A.combined),-log10(p.A_k.combined),fdr=fdr,method='median'),
                       apply(W.combined[,-(c(1:6, 14:20))],2,knockoff.threshold,fdr=fdr,offset=1),
                       # use MK.threshold() function to finalize corresponding thresholds
                       #MK.threshold(-log10(pValuesSTAAR.combined),-log10(pValuesSTAAR_k.combined),fdr=fdr,method='max'),
                       #MK.threshold(-log10(pValuesSTAAR.combined),-log10(pValuesSTAAR_k.combined),fdr=fdr,method='median'),
                       MK.threshold(-log10(pValuesACAT.combined),-log10(pValuesACAT_k.combined),fdr=fdr,method='max'),
                       MK.threshold(-log10(pValuesACAT.combined),-log10(pValuesACAT_k.combined),fdr=fdr,method='median'),
                       apply(W.combined[ , 16:18],2,knockoff.threshold,fdr=fdr,offset=1),
                       MK.threshold(-log10(pValuesACAT.combined),-log10(pValuesACAT_k_secondOrderComb.M),fdr=fdr,method='max'),
                       MK.threshold(-log10(pValuesACAT.combined),-log10(pValuesACAT_k_secondOrderComb.M),fdr=fdr,method='median'))
    Get.combineSelect <- function(x){W.combined[,x] >= thres.combine[x]}
    combineSelected <- sapply(1:ncol(W.combined),Get.combineSelect)
    combineSelected[is.na(combineSelected)] <- F
    # signal<-t(beta!=0)%*%window.matrix0!=0
    combineSignal <- c(t(beta!=0)%*%window.matrix0!=0, (beta!=0))
    
    combineSelected.BB <- as.vector(p.burden.combined<(0.05/ncol(X)))
    combineSelected.BD <- as.vector(p.dispersion.combined<(0.05/ncol(X)))
    combineSelected <- cbind(combineSelected, combineSelected.BB, combineSelected.BD)
    
    #BH FDR control
    combineSelected.BHB <- as.vector(p.adjust(p.burden.combined,method='fdr')<=fdr)
    combineSelected.BHD <- as.vector(p.adjust(p.dispersion.combined,method='fdr')<=fdr)
    combineSelected <- cbind(combineSelected, combineSelected.BHB, combineSelected.BHD)
    
    #selected<-cbind(selected,selected.RB,selected.RD)
    combinePower <- function(x){sum(combineSignal&x)/max(1,sum(combineSignal))}
    combineFdp <- function(x){sum((!combineSignal)&x)/max(1,sum(x))}
    
    combinedResult <- rbind(combinedResult,c(apply(combineSelected,2,combinePower),apply(combineSelected,2,combineFdp)))
  }
  if (replicateIndx==1) {combineResult.fdr <- combinedResult} else {combineResult.fdr <- combineResult.fdr + combinedResult}
  print(replicateIndx)
  combine.temp.print <- cbind(diseaseSnip, seq(0,0.2,by=0.01), combineResult.fdr/replicateIndx)
  print(combine.temp.print)
  
  combinedResult.detail <- rbind(combinedResult.detail, cbind(replicateIndx,h2,N.causal,seq(0,0.2,by=0.01),combinedResult))
  colnames(combinedResult.detail)<-c('replicate','h2','N.causal','fdr',
                                     'MSK_max_Burden_combine_power','MSK_max_Dispersion_combine_power','MSK_max_AllCauchy_combine_power',
                                     'MSK_median_Burden_combine_power','MSK_median_Dispersion_combine_power','MSK_median_AllCauchy_combine_power',
                                     'SK_Burden_combine_power','SK_Dispersion_combine_power','SK_AllCauchy_combine_power',
                                     'NormalK_Burden_combine_power','NormalK_Dispersion_combine_power','HMMK_Burden_combine_power','HMMK_Dispersion_combine_power',
                                     #'STAAR_max_combine_power', 'STAAR_median_combine_power',
                                     'STAAR_identity_max_combine_power', 'STAAR_identity_median_combine_power', 'STAAR_identity_single__combine_power',
                                     'NormalK_ACAT_combine_power', 'HMMK_ACAT_combine_power',
                                     'NormalK_multipleACAT_max_combine_power', 'NormalK_multipleACAT_median_combine_power',
                                     'Bonferroni_Burden_combine_power','Bonferroni_Dispersion_combine_power',
                                     'BH_Burden_combine_power','BH_Dispersion_combine_power', 
                                     'MSK_max_Burden_combine_fdp','MSK_max_Dispersion_combine_fdp','MSK_max_AllCauchy_combine_fdp',
                                     'MSK_median_Burden_combine_fdp','MSK_median_Dispersion_combine_fdp','MSK_median_AllCauchy_combine_fdp',
                                     'SK_Burden_combine_fdp','SK_Dispersion_combine_fdp','SK_AllCauchy_combine_fdp',
                                     'NormalK_Burden_combine_fdp','NormalK_Dispersion_combine_fdp','HMMK_Burden_combine_fdp','HMMK_Dispersion_combine_fdp',
                                     #'STAAR_max_combine_fdp', 'STAAR_median_combine_fdp',
                                     'STAAR_identity_max_combine_fdp', 'STAAR_identity_median_combine_fdp', 'STAAR_identity_single__combine_fdp',
                                     'NormalK_ACAT_combine_fdp', 'HMMK_ACAT_combine_fdp',
                                     'NormalK_multipleACAT_max_combine_fdp', 'NormalK_multipleACAT_median_combine_fdp',
                                     'Bonferroni_Burden_combine_fdp','Bonferroni_Dispersion_combine_fdp',
                                     'BH_Burden_combine_fdp','BH_Dispersion_combine_fdp')
  combinedResult.detail<-round(combinedResult.detail,digits=3)
  write.table(combinedResult.detail,paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/', location,'/indepCombineResult_knockoff_',out_type, '_',replicateIndx,'_',diseaseSnip,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}
