print('start')

args = commandArgs(TRUE)
print('job')
midResultFilePath = as.character(args[1])
#midResultFilePath="STATS_APOE4neg_Ages60plus_TOPMed_gen_source_UKB_ADSP_WGS_gen_source_Unrel_Rep_sDX"
print(c(midResultFilePath))

FDRTargetLevel = as.numeric(args[2])
#FDRTargetLevel=0.2
print(c(FDRTargetLevel))

library(data.table)
library(readr)
library(httr)
#library(liftOver)
library(plyr)
library(dplyr)
library(gdsfmt)
library(Matrix)
library(CMplot)
library(foreach)
library(doParallel)
source('/oak/stanford/groups/zihuai/UKB_KS/KnockoffScreen_AL_vMar16.R')
source('/oak/stanford/groups/zihuai/ESGWAS/ESGWAS.R')
source('/oak/stanford/groups/zihuai/ESGWAS_lasso/GhostKnockoff.R')
library(KnockoffScreen)

set.seed(123)


LD.dir<-paste0("/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/MetaMidResults/",midResultFilePath,"/")
LDmatrix.dir<-paste0("/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/LDMetaMidResults/",midResultFilePath,"/")
LD.names<-list.files(LD.dir)

################## GhostKnockoff filtering ################## 
resultsOverall<-NULL
for (i in 1:length(LD.names)) {
  print(i)
  gnomAD.data.all<-as.data.frame(fread(paste0(LD.dir,LD.names[i])))
  resultsOverall<-rbind.data.frame(resultsOverall,gnomAD.data.all[complete.cases(gnomAD.data.all),])
}
W_ghostScore<-MK.q.byStat(resultsOverall$kappa,resultsOverall$tau,5,Rej.Bound=10000)

#save mid stats for CM plot
out.dir<-paste0("/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/MidStats/",midResultFilePath,"/")
dir.create(file.path(out.dir), showWarnings=FALSE)

temp.filename<-paste0(out.dir,'Knockff_500_AD_results_Meta_singleLinkage','.csv')
midStats<-cbind.data.frame(resultsOverall,W_ghostScore)
write.table(midStats,temp.filename,col.names=T,row.names=F,sep='\t',quote=F)

ZScoreSelected<-W_ghostScore<=FDRTargetLevel
ZScoreSelected[is.na(ZScoreSelected)]<-F
SNPsRes<-cbind.data.frame(resultsOverall[ZScoreSelected,], W_ghostScore[ZScoreSelected])
colnames(SNPsRes)[ncol(SNPsRes)]<-"q"

#find neighboring SNPs from the selected blocks
#with correlation matrix
#when abs(cor)>corCutoff keep neighboring SNPs, otherwise remove; no filte: corCutoff=0
corCutoff <- 0.75
#extra filter with min Z-score from selected representative SNPs
minZScoreCutoff <- min(abs(SNPsRes$Z))

neighbourSNPs<-NULL
for (i in 1:length(LD.names)) {
  gnomAD.data.all<-as.data.frame(fread(paste0(LD.dir,LD.names[i])))
  if (length(unique(gnomAD.data.all$representhg38.chr))>1) {print(unique(gnomAD.data.all$representhg38.chr));print(i)}
  
  ifBlockExist<-which(SNPsRes$representhg38.chr %in% unique(gnomAD.data.all$representhg38.chr) & SNPsRes$representhg38.pos %in% unique(gnomAD.data.all$representhg38.pos))
  if (length(ifBlockExist)>0) {
    LDTemp<-as.matrix(fread(paste0(LDmatrix.dir,"matchedFilteredLD_",substring(LD.names[i], 40))),rownames=1)
    IDconsistency <- colnames(LDTemp) != paste0(gnomAD.data.all$hg19.pos,"_",gnomAD.data.all$ref,"_",gnomAD.data.all$alt)
    if (sum(IDconsistency)>0) {print(paste0(i," ID not consistent:", which(IDconsistency>0)))}
    
    for (j in ifBlockExist) {
      whichRepresentSNP<-which(gnomAD.data.all$hg38.chr==SNPsRes$representhg38.chr[j] & gnomAD.data.all$hg38.pos==SNPsRes$representhg38.pos[j])
      whichBlockExist<-which(gnomAD.data.all$representhg38.chr==SNPsRes$representhg38.chr[j] & gnomAD.data.all$representhg38.pos==SNPsRes$representhg38.pos[j])
      #find neighboring SNPs with abs cor>=cutoff
      corBlockTemp<-abs(LDTemp[whichBlockExist, whichRepresentSNP])
      filteredWhichBlockExist<-whichBlockExist[corBlockTemp>=corCutoff]
      
      qBlockTemp<-rep(NA,times=length(filteredWhichBlockExist))
      qBlockTemp[gnomAD.data.all$hg38.chr[filteredWhichBlockExist]==SNPsRes$representhg38.chr[j] &
                   gnomAD.data.all$hg38.pos[filteredWhichBlockExist]==SNPsRes$representhg38.pos[j]] <- SNPsRes[j,"q"]
      neighbourSNPs<-rbind.data.frame(neighbourSNPs,cbind.data.frame(gnomAD.data.all[filteredWhichBlockExist,],qBlockTemp,  corBlockTemp[corBlockTemp>=corCutoff]))
    }
  } 
}
#remove neighboring SNPs with min Z score
x1_sug <- as.data.frame(neighbourSNPs[abs(neighbourSNPs$Z) >= minZScoreCutoff,])
colnames(x1_sug)[ncol(x1_sug)] <- "withinClusterCorrelation"


################## Gene annotation ################## 
# find proximal genes
ref_seq<-read.table('/oak/stanford/groups/zihuai/GeneticsResources/refGene_hg38.txt')
colnames(ref_seq)[3]<-'chrom'
colnames(ref_seq)[5:6]<-c('txStart','txEnd')
colnames(ref_seq)[13]<-c('name2')

# find cS2G genes
topcS2GGene_allVariants <- read_csv("/oak/stanford/groups/zihuai/XinranQi/cS2G_UKBB/topcS2G_allVariants/topcS2GGene_allVariants.csv")
#locate cS2G genes for significant SNPs
precS2G_x1_sug <- topcS2GGene_allVariants[topcS2GGene_allVariants$hg19.pos>=min(x1_sug$hg19.pos) &
                                            topcS2GGene_allVariants$hg19.pos<=max(x1_sug$hg19.pos), ]
# find V2G genes
# Build query string
query_string = "
    query v2g($variantId: String!) {
    genesForVariant(variantId: $variantId) {
    gene {
        id
        symbol
        description
        bioType
        start 
        end
        chromosome
        }
    		functionalPredictions {
          aggregatedScore
          typeId
          sourceId
          tissues {
            maxEffectLabel
            maxEffectScore
            tissue {
              id
              name
            }
          }
        }
    		qtls {
          tissues {
            tissue {
              id
              name
            }
            quantile
            beta
            pval
          }
          typeId
          aggregatedScore
          sourceId
        }
        variant
        overallScore
        distances {
        sourceId
        aggregatedScore
        tissues {
            distance
          	tissue {
          	  id
          	} 
          	score 
          	quantile
        }
        }
        }
        }
"
# Set base URL of Genetics Portal GraphQL API endpoint
base_url <- "https://api.genetics.opentargets.org/graphql"

# fulfill gene annotations
x1_sug$proximalGene <- NA
x1_sug$cS2GGene <- x1_sug$cS2GMappingInfo <- x1_sug$cS2GScore <- NA
x1_sug$V2GGene <- x1_sug$V2GScore <- x1_sug$V2Geqtl <- x1_sug$V2Gpqtl <- x1_sug$V2Gsqtl <- NA
for (i in 1:nrow(x1_sug)) {
  # proximal
  ra <- x1_sug[i,]; colnames(ra)[1:2] <- c("hg38.chr","hg38.pos")
  ref_seq_t <- ref_seq[which(ref_seq$chrom==paste("chr",ra$hg38.chr,sep="")),]
  st_dist <- min(abs(ref_seq_t$txStart - ra$hg38.pos))
  st_dist_w <- which.min(abs(ref_seq_t$txStart - ra$hg38.pos))
  en_dist <- min(abs(ref_seq_t$txEnd - ra$hg38.pos))
  en_dist_w <- which.min(abs(ref_seq_t$txEnd - ra$hg38.pos))
  gene <- ""
  if (st_dist<=1000000 | en_dist<=1000000) {
    if (st_dist<=en_dist) {
      gene <- ref_seq_t[st_dist_w,"name2"]
    } else {
      gene <- ref_seq_t[en_dist_w,"name2"]
    }
  }
  x1_sug$proximalGene[i] <- ifelse(nchar(gene)==0,NA,as.character(gene))
  
  # cS2G
  cS2GPointer <- which(precS2G_x1_sug$hg19.chr==x1_sug$hg19.chr[i] & precS2G_x1_sug$hg19.pos==x1_sug$hg19.pos[i])
  if (length(cS2GPointer)>0) {
    cS2GTemp <- precS2G_x1_sug[cS2GPointer, c("topcS2GGene","cS2GScore","cS2GAnnotation")]; cS2GTemp <- cS2GTemp[!duplicated(cS2GTemp$topcS2GGene),]
    x1_sug$cS2GGene[i] <- cS2GTemp$topcS2GGene
    x1_sug$cS2GScore[i] <- cS2GTemp$cS2GScore
    x1_sug$cS2GMappingInfo[i] <- cS2GTemp$cS2GAnnotation
  }
  
  # V2G 
  tempAlleles<-x1_sug[i,]
  tempRefAlts<-if(nchar(tempAlleles$ref)>1 | nchar(tempAlleles$alt)>1) {
    c(paste0(tempAlleles$ref,"_",tempAlleles$alt), paste0(tempAlleles$alt,"_",tempAlleles$ref))
  } else {
    c(paste0(tempAlleles$ref,"_",tempAlleles$alt), paste0(tempAlleles$alt,"_",tempAlleles$ref),
      paste0(c("A","C","G","T")[tempAlleles$ref==c("T","G","C","A")],"_",c("A","C","G","T")[tempAlleles$alt==c("T","G","C","A")]),
      paste0(c("A","C","G","T")[tempAlleles$alt==c("T","G","C","A")],"_",c("A","C","G","T")[tempAlleles$ref==c("T","G","C","A")]))
  }
  hg38tempHits<-paste0(tempAlleles$hg38.chr,"_",tempAlleles$hg38.pos,"_",tempRefAlts)
  
  for (hitIndex in 1:length(hg38tempHits)) {
    variant_id <- hg38tempHits[hitIndex]#variant_id <- "1_55052794_A_G"
    # Set variables object of arguments to be passed to endpoint
    variables <- list("variantId"=variant_id)
    # Construct POST request body object with query string and variables
    post_body <- list(query=query_string, variables=variables)
    # Perform POST request
    r <- POST(url=base_url, body=post_body, encode='json')
    # Print first entry of L2G data to RStudio console
    V2GResult <- content(r)$data$genesForVariant
    if (length(V2GResult)>0) {break}
  }
  if(length(V2GResult)>0) {
    # choose V2G gene with highest overall score
    V2GResult <- V2GResult[[which.max(unlist(lapply(V2GResult, function(x){x$overallScore})))]]
    x1_sug$V2GGene[i] <- V2GResult$gene$symbol
    x1_sug$V2GScore[i] <- V2GResult$overallScore

    eqtlTemp <- sapply(V2GResult$qtls, function(x){x$sourceId=="eqtl"})
    x1_sug$V2Geqtl[i] <- ifelse(length(eqtlTemp)>0, ifelse(sum(eqtlTemp, na.rm=TRUE)>0, V2GResult$qtls[[which(eqtlTemp)]]$aggregatedScore,NA), NA)
    
    pqtlTemp <- sapply(V2GResult$qtls, function(x){x$sourceId=="pqtl"})
    x1_sug$V2Gpqtl[i] <- ifelse(length(pqtlTemp)>0, ifelse(sum(pqtlTemp, na.rm=TRUE)>0, V2GResult$qtls[[which(pqtlTemp)]]$aggregatedScore,NA), NA)
    
    sqtlTemp <- sapply(V2GResult$qtls, function(x){x$sourceId=="sqtl"})
    x1_sug$V2Gsqtl[i] <- ifelse(length(sqtlTemp)>0, ifelse(sum(sqtlTemp, na.rm=TRUE)>0, V2GResult$qtls[[which(sqtlTemp)]]$aggregatedScore,NA), NA)
  }
}


################## Single-cell RNA sequence differential expression ################## 
# scRNAseq gene level summary statistics
scRNAseq_geneLevel_cS2G_summaryStat <- read.delim("/oak/stanford/groups/zihuai/XinranQi/UKBGhost/scRNAseqGeneStatistics/cS2GGene/scRNAseq_geneLevel_cS2G_summaryStat.txt")
cS2GscRNAseq.summaryStat <- cbind.data.frame(scRNAseq_geneLevel_cS2G_summaryStat[,"cS2GGene"],
                                             apply(scRNAseq_geneLevel_cS2G_summaryStat[,c(5:18)], MARGIN=1, function(x){min(x, na.rm=TRUE)}),
                                             apply(scRNAseq_geneLevel_cS2G_summaryStat[,c(19:32)], MARGIN=1, function(x){min(x, na.rm=TRUE)}))
colnames(cS2GscRNAseq.summaryStat) <- c("gene","min_raw_pValue","min_adj_pValue")

scRNAseq_geneLevel_proximal_summaryStat <- read.delim("/oak/stanford/groups/zihuai/XinranQi/UKBGhost/scRNAseqGeneStatistics/proximalGene/scRNAseq_geneLevel_proximal_summaryStat.txt")
proximalscRNAseq.summaryStat <- cbind.data.frame(scRNAseq_geneLevel_proximal_summaryStat[,"proximalGene"],
                                                 apply(scRNAseq_geneLevel_proximal_summaryStat[,c(5:18)], MARGIN=1, function(x){min(x, na.rm=TRUE)}),
                                                 apply(scRNAseq_geneLevel_proximal_summaryStat[,c(19:32)], MARGIN=1, function(x){min(x, na.rm=TRUE)}))
colnames(proximalscRNAseq.summaryStat) <- c("gene","min_raw_pValue","min_adj_pValue")

scRNAseq.summaryStat <- rbind.data.frame(cS2GscRNAseq.summaryStat,proximalscRNAseq.summaryStat)
scRNAseq.summaryStat <- scRNAseq.summaryStat[!duplicated(scRNAseq.summaryStat$gene),]

x1_sug$cS2GscRNAseqMinP <- x1_sug$proximalscRNAseqMinP <- x1_sug$V2GscRNAseqMinP <- NA
# fulfillment of scRNAseq adjusted p-value
for (i in 1:nrow(x1_sug)) {
  x1_sug$cS2GscRNAseqMinP[i] <- ifelse(any(scRNAseq.summaryStat$gene==x1_sug$cS2GGene[i]),scRNAseq.summaryStat[scRNAseq.summaryStat$gene==x1_sug$cS2GGene[i],"min_adj_pValue"],NA)
  x1_sug$proximalscRNAseqMinP[i] <- ifelse(any(scRNAseq.summaryStat$gene==x1_sug$proximalGene[i]),scRNAseq.summaryStat[scRNAseq.summaryStat$gene==x1_sug$proximalGene[i],"min_adj_pValue"],NA)
  x1_sug$V2GscRNAseqMinP[i] <- ifelse(any(scRNAseq.summaryStat$gene==x1_sug$V2GGene[i]),scRNAseq.summaryStat[scRNAseq.summaryStat$gene==x1_sug$V2GGene[i],"min_adj_pValue"],NA)
}
x1_sug$W <- (x1_sug$kappa==0)*(x1_sug$tau)

x1_sug_output <- x1_sug[,c("hg38.chr","hg38.pos","hg19.chr","hg19.pos","ref","alt","AF","representhg38.chr","representhg38.pos", 
                           "qBlockTemp","W","pValue", "proximalGene","cS2GGene","cS2GMappingInfo","V2GGene","V2Gpqtl","V2Geqtl","V2Gsqtl", 
                           "proximalscRNAseqMinP","cS2GscRNAseqMinP","V2GscRNAseqMinP")]
colnames(x1_sug_output) <- c("hg38.chr","hg38.pos","hg19.chr","hg19.pos","REF","ALT","MAF","representhg38.chr","representhg38.pos", 
                             "q","W","p", "Proximal gene","cS2G >0.5 mapped gene","cS2G mapping info","Open Targets Genetics V2G gene","V2G aggregated pQTL score","V2G aggregated eQTL score","V2G aggregated sQTL score", 
                             "scRNAseq DEG minP of proximal gene","scRNAseq DEG minP of cS2G gene","scRNAseq DEG minP of V2G gene")
x1_sug_output <- apply(x1_sug_output,2,as.character)

#save mid stats for CM plot
tableOut.dir<-paste0("/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/FinalTables/VariantLevel/",midResultFilePath,"/")
dir.create(file.path(tableOut.dir), showWarnings=FALSE)
write.table(x1_sug_output,paste0(tableOut.dir,'stratifiedSummaryTable','.csv'),col.names=T,row.names=F,sep='\t',quote=F)


################## CM plot ################## 
# all variants
#x1 <- fread('/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/MidStats/STATS_APOE4neg_Ages60plus_TOPMed_gen_source_UKB_ADSP_WGS_gen_source_Unrel_Rep_sDX/Knockff_500_AD_results_Meta_singleLinkage.csv',header=T)
x1 <- midStats
x1 <- arrange(x1,hg38.chr,hg38.pos)
x1 <- as.data.frame(cbind(x1[,c("hg38.chr","hg38.pos","hg19.chr","hg19.pos","ref","alt","AF")],NA,x1[,c("Z","Z","kappa","tau","tau","W_ghostScore","pValue")]))
colnames(x1)<-c("CHR","BP","hg19.chr","hg19.pos","REF","ALT","MAF","isPermute","ES.Zscore_0","Schwartzentruber_etal","kappa","tau","W","q","_P")
x1$W<-(x1$kappa==0)*(x1$tau)
x1$SNP <- paste(x1$CHR,x1$BP,x1$REF,x1$ALT,sep=":")
x1 <- x1[which(!is.na(x1[,'CHR'])),]
#setorderv(x1, c("CHR","BP"))

# GhostKnockoff selected variants
x1_sug$KunklePValue <- x1_sug$MorenoPValue <- x1_sug$FinngenPValue <- NA
# for current CM plot, annotate each independent locus with proximal gene
x1_sug <- x1_sug[!is.na(x1_sug$proximalGene),]

x1_sug <- cbind(x1_sug[,c("hg38.chr","hg38.pos","hg19.chr","hg19.pos","ref","alt","AF")],NA,
                x1_sug[,c("Z","Z","kappa","tau","tau","qBlockTemp","pValue","representhg38.chr","representhg38.pos","cS2GGene","V2GGene","KunklePValue","MorenoPValue","FinngenPValue","proximalGene")])
colnames(x1_sug)<-c("CHR","BP","hg19.chr","hg19.pos","REF","ALT","MAF","isPermute","ES.Zscore_0","Schwartzentruber_etal","kappa","tau","W","q","_P","representhg38.chr","representhg38.pos","cS2G_GENE","V2G_GENE","KunklePValue","MorenoPValue","FinngenPValue","proximal_GENE")
x1_sug$W<-(x1_sug$kappa==0)*(x1_sug$tau)
#complete W values by replacing NA with representative SNPs' W values
x1_sug$Wcomplete<-x1_sug$W
x1_sug$qcomplete<-x1_sug$q
for (i in which(is.na(x1_sug$W))) {
  representSNP<-x1_sug[i,]
  representKappaTau<-x1[which(x1$CHR==representSNP$representhg38.chr & x1$BP==representSNP$representhg38.pos),c("kappa","tau","q")]
  x1_sug$Wcomplete[i]<-(representKappaTau$kappa==0)*(representKappaTau$tau)
  x1_sug$qcomplete[i]<-representKappaTau$q
}
x1_sug$SNP<-paste(x1_sug$CHR,x1_sug$BP,x1_sug$REF,x1_sug$ALT,sep=":")
#setorderv(x1_sug, c("CHR","BP"))
x1_sug <- arrange(x1_sug,CHR,BP)
#gene related to each selected SNPs, used for gene annotation of each independent locus
x1_sug$GENE<-x1_sug$proximal_GENE

# isolate SNPs and related data for selected SNPs
x1_sug$IND <- NA
x1_sug$TOP <- NA
x1_sug$RAst <- NA
x1_sug$RAen <- NA

loc_w <- 200000 #1000000
ra_w <- loc_w/2 #500000

# find independent loci
CHs <- unique(x1_sug$CHR)
inds <- 0
for (ch in CHs) {
  BPs <- x1_sug[which(x1_sug$CHR==ch),"BP"]
  BP_first <- BPs[1]
  for (bp in BPs) {
    if (bp==BP_first) {
      inds <- inds + 1
    } else {
      bp_sw <- bp - BP_prv
      if (bp_sw > loc_w) {
        inds <- inds + 1
      }
    }
    x1_sug[which(x1_sug$CHR==ch & x1_sug$BP==bp),"IND"] <- inds
    BP_prv <- bp
  }
}

# find top SNPs in independent loci
for (i in x1_sug$IND) {
  x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","SNP","W","_P","Wcomplete")]
  if (any(is.na(x1_sug_t$W))) {
    #first compare W values
    zeroPvalues <- if (sum(x1_sug_t$Wcomplete==max(x1_sug_t$Wcomplete))>1) {
      which(x1_sug_t$Wcomplete==max(x1_sug_t$Wcomplete))[which.min(x1_sug_t[x1_sug_t$Wcomplete==max(x1_sug_t$Wcomplete),"_P"])]
    } else {
      which.max(x1_sug_t$Wcomplete)
    }
    x1_sug_t[zeroPvalues,"TOP"] <- x1_sug_t[zeroPvalues,"SNP"]
    x1_sug[which(x1_sug$IND==i),c("TOP","SNP","W","_P","Wcomplete")] <- x1_sug_t
  } else {
    x1_sug_t[which.max(x1_sug_t$W),"TOP"] <- x1_sug_t[which.max(x1_sug_t$W),"SNP"]
    x1_sug[which(x1_sug$IND==i),c("TOP","SNP","W","_P","Wcomplete")] <- x1_sug_t
  }
}

# set range +/- extension for independent loci
for (i in x1_sug$IND) {
  x1_sug_t <- x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")]
  x1_sug_t[,"RAst"] <- x1_sug_t[1,"BP"] - ra_w
  x1_sug_t[,"RAen"] <- x1_sug_t[dim(x1_sug_t)[1],"BP"] + ra_w
  x1_sug[which(x1_sug$IND==i),c("TOP","BP","RAst","RAen")] <- x1_sug_t
}

#flag new loci
a <- unique(x1_sug[,c("representhg38.chr","representhg38.pos")])
#for each SNP, detect whether consistent with any of 3 independent studies
#if Consistent=TRUE, any of 3 p-values is significant
x1_sug$Consistent <- NA
for (i in 1:nrow(x1_sug)) {
  consistentTemp <- x1_sug[i,c("KunklePValue","MorenoPValue","FinngenPValue")]<=(0.05/nrow(a))
  consistentTemp <- consistentTemp[!is.na(consistentTemp)]
  x1_sug$Consistent[i] <- ifelse(length(consistentTemp)==0,NA,any(consistentTemp))
}

x1_sug_rank <- NULL
for (i in 1:nrow(a)) {
  x1_clusterTemp <- x1_sug[x1_sug$representhg38.chr==a$representhg38.chr[i] & x1_sug$representhg38.pos==a$representhg38.pos[i],]
  x1_clusterTemp <- cbind.data.frame(i,x1_clusterTemp[order(x1_clusterTemp[,"_P"], decreasing=FALSE),])
  x1_sug_rank <- rbind.data.frame(x1_sug_rank,x1_clusterTemp)
}
colnames(x1_sug_rank)[1] <- "clusterIndex"

#for each gene, detect whether new discovery
geneNewDiscovery <- rep(NA, times=length(unique(x1_sug$GENE)))
names(geneNewDiscovery) <- unique(x1_sug$GENE)
for (i in 1:length(unique(x1_sug$GENE))) {
  newHitTemp <- x1_sug_rank$Consistent[x1_sug_rank$GENE==unique(x1_sug$GENE)[i]]
  names(newHitTemp) <- x1_sug_rank$clusterIndex[x1_sug_rank$GENE==unique(x1_sug$GENE)[i]]
  newHitTemp <- newHitTemp[!is.na(newHitTemp)]
  geneNewDiscovery[i] <- ifelse(length(newHitTemp)==0,NA,any(newHitTemp))
}

#for a gene linked to variant, if any cluster level consistency with any of 3 independent studies, then it NOT new_hit
#NEW_hit=FALSE meaning consistent with any 3 studies at cluster level
#NEW_hit=TRUE meaning not consistent with any 3 studies at cluster level or NA
x1_sug$NEW_hit <- NA
for(i in 1:nrow(x1_sug)){
  newHitTemp <- geneNewDiscovery[names(geneNewDiscovery)==x1_sug$GENE[i]]
  x1_sug$NEW_hit[i] <- ifelse(is.na(newHitTemp),"Y",!newHitTemp)
}

# extract signal from all variants for respective independent loci (for plotting)
#x1_sig <- x1_sug[c(which(x1_sug[,'q']<=0.05),which(is.na(x1_sug$q))[x1_sug[is.na(x1_sug$q),"_P"]<1e-10]),]
x1_sig <- x1_sug[x1_sug[,'qcomplete']<=FDRTargetLevel, ]
#setorderv(x1_sig, c("CHR","BP"))
x1_sig <- arrange(x1_sig,CHR,BP)

signal <- NULL
for (t in x1_sig[which(!is.na(x1_sig$TOP)),"TOP"]) {
  ra <- x1_sig[which(x1_sig$TOP==t),c("CHR","RAst","RAen","GENE","q","W","IND","SNP","Wcomplete","qcomplete","cS2G_GENE")]
  xs <- x1_sug[which(x1_sug$CHR==ra$CHR & (x1_sug$BP>=ra$RAst) & (x1_sug$BP<=ra$RAen)), c("SNP","q","W","BP","CHR","NEW_hit","GENE","_P","Wcomplete","cS2G_GENE")]
  if (ra$qcomplete<=FDRTargetLevel) {xs$col <- "purple"}
  geneCounts <- table(xs$GENE)
  #keep first 3 cS2G genes
  numOfMaxDuplicateGene <- ifelse(length(names(geneCounts))>=3,3,length(names(geneCounts)))
  indexOfMaxDuplicateGene <- order(geneCounts, decreasing=TRUE)[1:numOfMaxDuplicateGene]
  maxDuplicates <- geneCounts[indexOfMaxDuplicateGene]
  maxDuplicateName <- names(geneCounts)[indexOfMaxDuplicateGene]
  
  if (length(unique(xs$GENE))!=nrow(xs)) {
    #for gene with more than 1 times of occurance, choose one with biggest W
    keepRowIndex <- NULL
    for (geneIndx in names(geneCounts)[geneCounts>1]) {
      xsTemp <- xs[which(xs$GENE==geneIndx),]
      if (any(is.na(xsTemp$W))) {
        #first compare W values
        zeroPvalues <- if (sum(xsTemp$Wcomplete==max(xsTemp$Wcomplete))>1) {
          which(xsTemp$Wcomplete==max(xsTemp$Wcomplete))[which.min(xsTemp[xsTemp$Wcomplete==max(xsTemp$Wcomplete),"_P"])]
        } else {
          which.max(xsTemp$Wcomplete)
        }
      } else {
        zeroPvalues <- which.max(xsTemp$W)
      }
      keepRowIndex <- c(keepRowIndex, which(xs$GENE==geneIndx)[zeroPvalues])
    }
    xsAllGene <- xs[c(keepRowIndex,which(xs$GENE %in% names(geneCounts)[geneCounts==1])),]
    xs <- xsAllGene[which.min(xsAllGene[,"_P"]),]
    xsGeneCount <- geneCounts[names(geneCounts)==xs$GENE]
  } else {
    xs <- xs[which.min(xs[,"_P"]),]
    xsGeneCount <- 1
  }
  #xs$text_col <- sapply(1:nrow(xs), function(x){ifelse(is.na(xs[x,"cS2G_GENE"]),"#009E73","black")})
  xs$text_col <- sapply(1:nrow(xs), function(x){ifelse(xs[x,"NEW_hit"]=="Y","red","black")})
  xs$pch <- 19
  xs$text <- xs$GENE # paste0(xs$GENE,"(",xsGeneCount,")")
  xs$text1 <- maxDuplicateName[1] # paste0(maxDuplicateName[1],"(",maxDuplicates[1],")")
  signal <- rbind.data.frame(signal,xs)
}
signal$text_col<-"black"

# keep signal with highlights
# complete W with same cluster representative SNP's W
x1Supplement <- NULL
for (i in which(!complete.cases(signal))) {
  representSNP <- x1_sug[x1_sug$SNP==signal$SNP[i],]
  signal[i,c("q","W")] <- x1[which(x1$CHR==representSNP$representhg38.chr & x1$BP==representSNP$representhg38.pos),c("q","W")]
  x1SupplementTemp <- x1[which(x1$CHR==representSNP$representhg38.chr & x1$BP==representSNP$representhg38.pos),]
  x1SupplementTemp[,c("CHR","BP","hg19.chr","hg19.pos","REF","ALT","MAF","ES.Zscore_0","Schwartzentruber_etal","_P","SNP")] <- representSNP[,c("CHR","BP","hg19.chr","hg19.pos","REF","ALT","MAF","ES.Zscore_0","Schwartzentruber_etal","_P","SNP")]
  x1Supplement <- rbind.data.frame(x1Supplement,x1SupplementTemp)
}
signal_topp <- signal

x1 <- rbind.data.frame(x1,x1Supplement)
#setorderv(x1, c("CHR","BP"))
x1 <- arrange(x1,CHR,BP)

# thresholds
T010 <- min(x1[x1$q<=FDRTargetLevel,"W"])#; T005 <- x1_sug[which.min(abs(x1_sug$q-0.05)),"W"]
ths <- c(T010)
#ths <- c(T010,T005)
# data points in CM plot
x1t <- x1[,c("SNP","CHR","BP","W")] 
x1t[which(x1t$W>100),"W"] <- 100
main.text<-midResultFilePath

CMPlotOut_dir <- paste0("/oak/stanford/groups/zihuai/XinranQi/stratifiedGhostKnockoff/CMPlots/",midResultFilePath,"/")
dir.create(file.path(CMPlotOut_dir), showWarnings=FALSE)
setwd(CMPlotOut_dir)

CMplot(x1t, plot.type="m", LOG10=FALSE, col=c("grey30","grey60"), ylab="W statistic",bin.range=c(0,1000),
       chr.den.col=c("darkgreen", "yellow", "red"), main=main.text,
       threshold=ths, threshold.lty=c(2,2), threshold.lwd=c(1,1), threshold.col=c("black","red"),
       highlight=signal_topp$SNP, highlight.cex=1, 
       highlight.col=signal_topp$col, highlight.text.col=signal_topp$text_col, highlight.text=signal_topp$text,
       signal.col=c("cornflowerblue"),signal.cex=c(1),main.cex=1,
       file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)
