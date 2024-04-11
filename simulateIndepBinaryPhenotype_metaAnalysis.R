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
#numIndependentDiseased= numIndependentControl=25000 
#numOfSelectedSnips=1000 
#diseaseSnip=10 
#location="SchemeA/D/meta" 
#h2=2.5; N.causal=10; out_type="D";theta=NA

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
library(doMC)
library(irlba)
library(ff)
library(pROC)
library(aod)
library(lmtest)
library(poolr)

source('/oak/stanford/groups/zihuai/XinranQi/Project1/simulate10PedigreeRecombine0.R')
source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
library(KnockoffScreen)

data("SKAT.haplotypes")
SKATHaplotypes <- SKAT.haplotypes$Haplotype[ , SKAT.haplotypes$SNPInfo$FREQ1>0.001]
positionChromosome <- SKAT.haplotypes$SNPInfo$CHROM_POS[SKAT.haplotypes$SNPInfo$FREQ1>0.001]

set.seed(123)

singleResult.detail <- c()
for (replicateIndx in 1:1000){
  print(replicateIndx)
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
  caseIndx <- sample(which(y==1), 1000, replace=FALSE)
  controlIndx <- sample(which(y==0), 9000, replace=FALSE)
  G <- G[c(caseIndx, controlIndx), ]
  X1 <- X1[c(caseIndx, controlIndx)]
  y <- y[c(caseIndx, controlIndx)]
  
  #summary stats -- always use shrinkage
  cor.G <- matrix(as.numeric(corpcor::cor.shrink(G,verbose=F)), nrow=ncol(G))
  # fit null model
  fit.prelim.Ghost <- GhostKnockoff.prelim(cor.G,M=5,method='sdp',corr_max=0.75)

  # apply other alternative methods
  result.prelim <- KS.prelim(y,X=X1,out_type=out_type)
  mu<-result.prelim$nullglm$fitted.values;Y.res<-result.prelim$Y-mu
  p.add<-as.matrix(Get.p(G,result.prelim))

  # Ghostknockoff score test
  n<-length(Y.res)
  Zscore_0<-qnorm(p.add/2,lower.tail=F)*sign(cov(as.matrix(G),Y.res)) #convert p-values to z-scores
  GK.stat<-GhostKnockoff.fit(Zscore_0=Zscore_0,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  GK.filter <- MK.statistic(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
  W_ghostScore <- MK.q.byStat(GK.filter[,'kappa'],GK.filter[,'tau'],5,Rej.Bound=10000)

  #meta p-values
  #half y==1 & half y==0
  oneSampleIndex<-1:(n/2)
  otherSampleIndex<-(n/2+1):n
  #oneSampleIndex<-sort(sample(n, n/2, replace=FALSE))
  #otherSampleIndex<-setdiff(1:n,oneSampleIndex)
  
  result.prelim1<-KS.prelim(y[oneSampleIndex],X=X1[oneSampleIndex],out_type=out_type)
  mu1<-result.prelim1$nullglm$fitted.values;Y.res1<-result.prelim1$Y-mu1
  p.addMeta1<-as.matrix(Get.p(G[oneSampleIndex,],result.prelim1)) #score test
  
  result.prelim2<-KS.prelim(y[otherSampleIndex],X=X1[otherSampleIndex],out_type=out_type)
  mu2<-result.prelim2$nullglm$fitted.values;Y.res2<-result.prelim2$Y-mu2
  p.addMeta2<-as.matrix(Get.p(G[otherSampleIndex,],result.prelim2)) #score test
  
  p.addMeta<-sapply(1:length(p.add), function(x){fisher(c(p.addMeta1[x],p.addMeta2[x]), adjust="none")$p})
  
  Zscore_0Meta<-qnorm(p.addMeta/2,lower.tail=F)*sign(cov(as.matrix(G),rbind(Y.res1,Y.res2))) #convert p-values to z-scores
  GK.statMeta<-GhostKnockoff.fit(Zscore_0=Zscore_0Meta,N.effect=n,fit.prelim=fit.prelim.Ghost,method='marginal')
  GK.filterMeta <- MK.statistic(GK.statMeta$T_0[[1]],GK.statMeta$T_k[[1]])
  W_ghostScoreMeta <- MK.q.byStat(GK.filterMeta[,'kappa'],GK.filterMeta[,'tau'],5,Rej.Bound=10000)
  
  # knockoff-based inference for Z scores with/out sample relatedness adjustment
  W.ZScores <- cbind( W_ghostScore, W_ghostScoreMeta); W.ZScores[is.na(W.ZScores)]<-0
  
  # single variable knockoff-based inference methods
  singleResult <- c()
  for (fdr in seq(0,0.2,by=0.01)) {
    Get.ZScoreSelect <- function(x){W.ZScores[,x] <= fdr}
    ZScoreSelected <- sapply(1:ncol(W.ZScores), Get.ZScoreSelect)
    ZScoreSelected[is.na(ZScoreSelected)] <- F

    singleSignal <- (beta!=0)
    
    singlePower <- function(x) {sum(singleSignal&x)/max(1,sum(singleSignal))}
    singleFdp <- function(x) {sum((!singleSignal)&x)/max(1,sum(x))}
    
    singleResult <- rbind(singleResult, c(apply(ZScoreSelected,2,singlePower), apply(ZScoreSelected,2,singleFdp)))
  }
  if (replicateIndx==1) {singleResult.fdr <- singleResult} else {singleResult.fdr <- singleResult.fdr + singleResult}
  print(replicateIndx)
  single.temp.print <- cbind(N.causal, seq(0,0.2,by=0.01), singleResult.fdr/replicateIndx)
  print(single.temp.print)
  
  singleResult.detail<-rbind(singleResult.detail, cbind(replicateIndx,h2,N.causal,seq(0,0.2,by=0.01),singleResult))
  colnames(singleResult.detail)<-c('replicate','h2','N.causal','fdr',
                                   'Ghostknockoff_score_power', 'Ghostknockoff_metaScore_power', 
                                   'Ghostknockoff_score_fdp', 'Ghostknockoff_metaScore_fdp')
  singleResult.detail <- round(singleResult.detail,digits=5)
  write.table(singleResult.detail,paste0('/oak/stanford/groups/zihuai/XinranQi/Project1/', location,'/IndepSingleResult_knockoff_',out_type, '_',replicateIndx,'_',theta,'_detail.txt'),row.names=F,col.names=T,quote=F,sep='\t')
}
