print('start')

args = commandArgs(TRUE)
print('job')
job_index = as.numeric(args[1])
print(c(job_index))

job_seq = as.numeric(args[2])
print(c(job_seq,job_index))

#ml R/4.0.2
library(liftOver)
library(data.table)
library(gdsfmt)
library(Matrix)
source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
library(KnockoffScreen)

set.seed(123)

#LD.dir<-'/oak/stanford/groups/zihuai/gnomAD/LD_Scores/LD_matrices/'
LD.dir<-'/oak/stanford/groups/zihuai/gnomAD/LD_Scores/nearly_independent_Beriza/complete_LD_matrices_csv/'
LD.names<-list.files('/oak/stanford/groups/zihuai/gnomAD/LD_Scores/nearly_independent_Beriza/complete_LD_matrices_csv/')
LD.names<-LD.names[LD.names!="partition_summary.csv"]
out.dir<-'/oak/stanford/groups/zihuai/XinranQi/UKBGhost/MetaMidResults/'

Zscores<-read.delim("/oak/stanford/groups/zihuai/ESGWAS/GWAS_SummaryStat/Belloy/GCST90012878_Schwartzentruber_2021_UKB_GRCh38_MAF_0.01_Zscore.txt")
#Zscores <- read.delim("/oak/stanford/groups/zihuai/ESGWAS/GWAS_SummaryStat/Belloy/Belloy_2021_June_GWAS_AD_UKB_EUR_GRCh38_N438080_LMM_scDX_MAF_0.01_Zscore.txt")
#Zscores <- read.delim("/oak/stanford/groups/zihuai/ESGWAS/GWAS_SummaryStat/Belloy/Belloy_2021_June_GWAS_AD_UKB_EUR_GRCh38_N438080_LMM_uDX_MAF_0.01_Zscore.txt")
Zscores$pValue<-sapply(Zscores$Z, function(x){2*min(pnorm(q=x,lower.tail=FALSE),pnorm(q=x,lower.tail=TRUE))})
 
dictionaryAllSNP <- fread("/oak/stanford/groups/zihuai/XinranQi/cS2G_V2G_results/Results/MidStats/gnomadQC_remove.txt")
#the biggest file exceeds memory limit, thus divide it only with object function: maximize (L-x)^2+x^2
#ffinf<-file.info(dir(LD.dir, full.names=TRUE), extra_cols=FALSE)
#utils:::format.object_size(ffinf$size[which.max(ffinf$size)], "auto")
#which(ffinf$size==sort(ffinf$size,decreasing=TRUE)[2])

#library(readr)
#GhostHits<-read_csv("/oak/stanford/groups/zihuai/XinranQi/GhostKnockoffSussie.stats.sug_hits.info.csv")
#GhostHitsAPOE<-GhostHits[GhostHits$CHR==19,]
#ZscoresAPOE<-Zscores[Zscores$BP%in%GhostHitsAPOE$BP & Zscores$CHR==19,]
#cat(sum(ZscoresAPOE$pValue<5e-8)," out of ", nrow(ZscoresAPOE), " APOE SNPs have p-values smaller than 5e-8.")
#write_csv(ZscoresAPOE,'/oak/stanford/groups/zihuai/XinranQi/UKBGhost/ZScoresAPOE.csv')

job_seq.tag<-c(seq(1,length(LD.names),by=100),length(LD.names))
LD.names<-LD.names[job_seq.tag[job_seq]:job_seq.tag[job_seq+1]]
N_set<-length(LD.names)
options(scipen=999)
set<-job_index

#job_seq=7;job_index=2

if(set>N_set){print ("analysis already done")}else{
  print(c(job_seq,job_index))
  
  gnomAD.data.all<-as.data.frame(fread(paste0(LD.dir,LD.names[job_index])))
  variant.info<-as.data.frame(gnomAD.data.all[,2:6])
  chr<-variant.info[,'chr']
  pos<-variant.info[,'pos']

  #hg19 coordinates for LD matrix
  df<-cbind(data.frame(paste0('chr',chr)),pos,pos)
  colnames(df)<-c('chr','start','end')
  temp.Granges<-makeGRangesFromDataFrame(df)
  #chain <- import.chain("/oak/stanford/groups/zihuai/GeneticsResources/LiftOver/hg38ToHg19.over.chain")
  chain <- import.chain("/oak/stanford/groups/zihuai/GeneticsResources/LiftOver/hg19ToHg38.over.chain")
  converted.df<-data.frame(liftOver(temp.Granges,chain))
  converted.chr<-as.numeric(sub('chr','',converted.df[match(1:length(pos),converted.df[,1]),'seqnames']))
  converted.pos<-converted.df[match(1:length(pos),converted.df[,1]),'start']
  
  info<-cbind(converted.chr,converted.pos,variant.info)
  colnames(info)[1:4]<-c('hg38.chr','hg38.pos','hg19.chr','hg19.pos')
  
  #load Zscores
  print(paste0("LD matrix of chromosome ", unique(info$hg38.chr)))
  uniqueCHR<-if(any(is.na(unique(info$hg38.chr)))) {unique(info$hg38.chr)[!is.na(unique(info$hg38.chr))]} else {unique(info$hg38.chr)}
  
  if(length(uniqueCHR)==1) {
    Zscores.chr<-Zscores[Zscores[,'CHR']==uniqueCHR,]
    temp.Zscores<-Zscores.chr[which(Zscores.chr[,'BP'] >= min(info$hg38.pos,na.rm=TRUE) & Zscores.chr[,'BP'] <= max(info$hg38.pos,na.rm=TRUE)),]
    
    include.index<-match(intersect(temp.Zscores[,'BP'],info$hg38.pos),info$hg38.pos)
    gnomAD.data.all<-gnomAD.data.all[include.index,c(2:6,(include.index+(ncol(gnomAD.data.all)-nrow(gnomAD.data.all))))]
    infoMatch<-info[include.index,]
    
    include.index<-match(intersect(temp.Zscores[,'BP'],info$hg38.pos),temp.Zscores[,'BP'])
    temp.Zscores<-temp.Zscores[include.index,]
  } else {
    Zscores.chr<-Zscores[Zscores[,'CHR'] %in% uniqueCHR,]
    temp.Zscores<-Zscores.chr[which(Zscores.chr[,'BP'] >= min(info$hg38.pos,na.rm=TRUE) & Zscores.chr[,'BP'] <= max(info$hg38.pos,na.rm=TRUE)),]
    
    include.index<-match(intersect(paste0(temp.Zscores[,'CHR'],"_",temp.Zscores[,'BP']),paste0(info$hg38.chr,"_",info$hg38.pos)),paste0(info$hg38.chr,"_",info$hg38.pos))
    gnomAD.data.all<-gnomAD.data.all[include.index,c(2:6,(include.index+(ncol(gnomAD.data.all)-nrow(gnomAD.data.all))))]
    infoMatch<-info[include.index,]
    
    include.index<-match(intersect(paste0(temp.Zscores[,'CHR'],"_",temp.Zscores[,'BP']),paste0(info$hg38.chr,"_",info$hg38.pos)),paste0(temp.Zscores[,'CHR'],"_",temp.Zscores[,'BP']))
    temp.Zscores<-temp.Zscores[include.index,]
  }
  
  #match reference & alternative alleles
  alleleMatch<-signZscores<-numeric(nrow(temp.Zscores))
  for (matchIndx in 1:nrow(temp.Zscores)) {
    tempAlleles<-gnomAD.data.all[infoMatch$hg38.chr==temp.Zscores[matchIndx,"CHR"] & infoMatch$hg38.pos==temp.Zscores[matchIndx,"BP"], c("ref","alt")]
    tempAlleles<-if(nchar(tempAlleles$ref)>1 | nchar(tempAlleles$alt)>1) {
      rbind(tempAlleles, unlist(tempAlleles)[c(2,1)])
    } else {
      rbind(tempAlleles, unlist(tempAlleles)[c(2,1)],
            c(c("A","C","G","T")[tempAlleles$ref==c("T","G","C","A")], c("A","C","G","T")[tempAlleles$alt==c("T","G","C","A")]),
            c(c("A","C","G","T")[tempAlleles$alt==c("T","G","C","A")], c("A","C","G","T")[tempAlleles$ref==c("T","G","C","A")]))
    }
    tempAllelesWithin<-apply(tempAlleles, 1, function(x){x[1]==temp.Zscores[matchIndx,"REF"] & x[2]==temp.Zscores[matchIndx,"ALT"]})
    names(tempAllelesWithin)<-NULL
    alleleMatch[matchIndx]<-any(tempAllelesWithin)
    signZscores[matchIndx]<-ifelse(any(tempAllelesWithin),ifelse(min(which(tempAllelesWithin))%in%c(2,4),-1,1),NA)
  }
  temp.Zscores<-cbind.data.frame(infoMatch[,c('hg38.chr','hg38.pos','hg19.chr','hg19.pos','ref','alt','AF')],temp.Zscores$Z*signZscores,temp.Zscores$pValue)
  colnames(temp.Zscores)[c(ncol(temp.Zscores)-1,ncol(temp.Zscores))]<-c("Z","pValue")
  temp.Zscores<-temp.Zscores[which(alleleMatch>0),]
  
  gnomAD.data.all<-gnomAD.data.all[which(alleleMatch>0),c(1:5,(which(alleleMatch>0)+ncol(gnomAD.data.all)-nrow(gnomAD.data.all)))]
 
  # gnomAD filtering to filter out variants that are not passing gnomAD QC
  removeIndx<-rep(NA,times=nrow(temp.Zscores))
  for (removei in 1:nrow(temp.Zscores)) {
    print(removei)
    removeIndx[removei]<-any(paste0("chr",temp.Zscores$hg19.chr[removei])==dictionaryAllSNP$CHR & temp.Zscores$hg19.pos[removei]==dictionaryAllSNP$BP)
  }
  temp.Zscores<-temp.Zscores[which(removeIndx==FALSE),]
  gnomAD.data.all<-gnomAD.data.all[which(removeIndx==FALSE),c(1:5,(which(removeIndx==FALSE)+ncol(gnomAD.data.all)-nrow(gnomAD.data.all)))]
   
  #data frame indicating cluster using representative SNPs
  clusterRepresent<-temp.Zscores[,c('hg38.chr','hg38.pos','hg19.chr','hg19.pos','ref','alt','AF','Z','pValue')]
  clusterRepresent<-cbind.data.frame(clusterRepresent,NA,NA,NA,NA)
  colnames(clusterRepresent)[(ncol(clusterRepresent)-3):ncol(clusterRepresent)]<-c('representhg38.chr','representhg38.pos',"kappa","tau")

  #correlation matrix
  cor.G<-gnomAD.data.all[,-(1:5)]
  cor.G<-as.matrix(cor.G)
  cor.G<-cor.G+t(cor.G)
  diag(cor.G)<-diag(cor.G)/2
  cor.G<-as.matrix(forceSymmetric(cor.G))
  
  #data adaptive tuning for hyper-parameter "max.size"
  #withinCorLess<-FALSE
  #maxSizeRealization<-500
  #while (!withinCorLess) {
  #  block.fit<-divide.sdp(cor.G,max.size=maxSizeRealization)
  #  block.num<-length(table(block.fit$clusters))
  #  quantile75s<-rep(NA,block.num)
  #  for(k in 1:block.num){
      #across cluster
  #    across.block.cor<-as.matrix(cor.G[which(block.fit$clusters==k),which(block.fit$clusters!=k)])
  #    acrossEntries<-abs(as.vector(across.block.cor))
  #    quantile75s[k]<-quantile(acrossEntries,probs=0.75,na.rm=TRUE)
  #  }
  #  withinCorLess<-mean(quantile75s,na.rm=TRUE)<0.05
  #  maxSizeRealization<-ifelse(withinCorLess,maxSizeRealization,maxSizeRealization+500)
    #compare with maximum # of variants
  #  if(maxSizeRealization>(nrow(cor.G)/2)) {withinCorLess<-TRUE; maxSizeRealization<-round(nrow(cor.G)/2,digits=0)}
    #print(c(mean(quantile75s,na.rm=TRUE),maxSizeRealization,nrow(cor.G)))
  #}
  block.fit<-divide.sdp(cor.G,max.size=500)
  block.num<-length(table(block.fit$clusters))
  
  for(k in 1:block.num){
    print(k)
    block.index<-which(block.fit$clusters==k)
    gnomAD.data<-gnomAD.data.all[block.index,c(1:5,block.index+5)]
    block.cor.G<-gnomAD.data[,-(1:5)]
    block.cor.G<-as.matrix(block.cor.G)
    block.cor.G<-block.cor.G+t(block.cor.G)
    diag(block.cor.G)<-diag(block.cor.G)/2
    block.cor.G<-as.matrix(forceSymmetric(block.cor.G))

    eigen.fit<-eigen(block.cor.G)
    newEig<-ifelse(eigen.fit$values<1e-5, 1e-5, eigen.fit$values)
    newMat<-eigen.fit$vectors %*% (newEig*t(eigen.fit$vectors))
    # normalize modified matrix eqn 6 from Brissette et al 2007
    newMat<-newMat/sqrt(diag(newMat) %*% t(diag(newMat)))
    block.cor.G<-newMat
    
    Sigma<-block.cor.G
    corr_max<-0.75
    Sigma.distance<-as.dist(1-abs(Sigma))
    if(ncol(Sigma)>1){
      fit<-hclust(Sigma.distance, method="single")
      clusters<-cutree(fit, h=1-corr_max)
    }else{clusters<-1}
    cluster.index<-c()
    for(j in 1:max(clusters)){
      temp.cor<-block.cor.G[clusters==j,clusters==j,drop=F]
      cluster.index<-c(cluster.index,which(clusters==j)[which.max(colSums(abs(temp.cor)))])
    }
    
    #record cluster with representative SNPs & kappa & tau & hg38/19 CHR BP
    representativeSNPIndx<-numeric(length(clusters))
    for (clusterIndx in sort(unique(clusters))) {
      representativeSNPIndx[clusters==clusters[cluster.index[clusterIndx]]]<-cluster.index[clusterIndx]
    }
    clusterRepresent[which(block.fit$clusters==k),(ncol(clusterRepresent)-3):(ncol(clusterRepresent)-2)]<-temp.Zscores[which(block.fit$clusters==k)[representativeSNPIndx],c("hg38.chr","hg38.pos")]
    
    gnomAD.data<-gnomAD.data[cluster.index,c(1:5,cluster.index+5)]
    block.cor.G<-block.cor.G[cluster.index,cluster.index,drop=F]

    block.MAF<-as.matrix(gnomAD.data[,5])
    block.MAF[block.MAF>0.5]<-1-block.MAF[block.MAF>0.5]
    block.REF<-as.matrix(gnomAD.data[,'ref'])
    block.ALT<-as.matrix(gnomAD.data[,'alt'])
    block.id<-paste0(as.matrix(gnomAD.data[,'chr']),':',as.matrix(gnomAD.data[,'pos']),'-',as.matrix(gnomAD.data[,'ref']),'-',as.matrix(gnomAD.data[,'alt']))
    colnames(block.cor.G)<-block.id

    block.cor.G<-as.matrix(forceSymmetric(block.cor.G))
    fit.prelim.Ghost<-GhostKnockoff.prelim(block.cor.G,M=5,method='sdp',corr_max=0.75)
    
    #remaining Z scores corresponding to cluster.index
    tempZscoresCluster<-temp.Zscores[which(block.fit$clusters==k)[cluster.index],,drop=F]
    tempZscoresCluster$Z[is.na(tempZscoresCluster$Z)]<-0
    
    # Ghostknockoff from BOLT-LMM p values
    GK.stat<-GhostKnockoff.fit(Zscore_0=tempZscoresCluster$Z,N.effect=nrow(tempZscoresCluster),fit.prelim=fit.prelim.Ghost,method='marginal')
    GK.filter <- MK.statistic(GK.stat$T_0[[1]],GK.stat$T_k[[1]])
    #record kappa & tau
    clusterRepresent[which(block.fit$clusters==k)[cluster.index],(ncol(clusterRepresent)-1):ncol(clusterRepresent)]<-cbind(GK.filter[,'kappa'],GK.filter[,'tau'])
  }
  temp.filename<-paste0(out.dir,'meta_mid_statistics_500_',sub('.csv.*','',LD.names[job_index]),'.csv')
  write.table(clusterRepresent,temp.filename,col.names=T,row.names=F,sep='\t',quote=F)
}
