library(dplyr)
library(ggplot2)
library("ggpubr")
library(patchwork)
library(stringr)
library(ggrepel)
library(gridExtra)
library(ggthemes)
library(grid)

############# related H0 ################
GIF <- function(PValues) {
  # For p-values, calculate chi-squared statistic
  chisq <- qchisq(PValues,1)
  return(qchisq(0.5,1)/median(chisq))
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

aggregatePValue <- function(thetaValue, file_path) {
  files <- list.files(file_path)
  files <- files[sapply(1:length(files), function(x){substrRight(files[x],4)==".txt"})]
  
  pValueTrace <- NULL
  for (i in 1:length(files)) {
    filei <- read.delim(paste0(file_path, "/", files[i]))
    pValueTrace <- rbind(pValueTrace, filei)
  }
  pValueTrace <- cbind(pValueTrace, thetaValue)
  # output
  return(as.data.frame(pValueTrace))
}

ACATSTAAR_0.4_0_2_detail <- aggregatePValue(0, "/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=0")
ACATSTAAR_0.4_4_2_detail <- aggregatePValue(4, "/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=4")
ACATSTAAR_0.4_7.5_2_detail <- aggregatePValue(7, "/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=7")
ParameterEst <- rbind.data.frame(ACATSTAAR_0.4_0_2_detail, ACATSTAAR_0.4_4_2_detail, ACATSTAAR_0.4_7.5_2_detail)

ACATSTAAR_0.4_0_2_detail <- aggregatePValue(0, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeA/D/theta=0")
ACATSTAAR_0.4_4_2_detail <- aggregatePValue(4, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeA/D/theta=4")
ACATSTAAR_0.4_7_2_detail <- aggregatePValue(7, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeA/D/theta=7")
ParameterEst <- rbind.data.frame(ACATSTAAR_0.4_0_2_detail,  ACATSTAAR_0.4_4_2_detail, ACATSTAAR_0.4_7_2_detail)

colnames(ParameterEst) <- c("window", "single", "Class", "theta")
ParameterEst <- ParameterEst %>% arrange(Class)


parameterEstimateComparison <- function(Thetas, VarianceValue) {
  windowPlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubset <- as.data.frame(ParameterEst[ParameterEst$theta==Thetas[thetaIndx], c(1,3)])
    dataSubset <- dataSubset[complete.cases(dataSubset), ]
    dataSubset$window <- as.numeric(dataSubset$window)
    dataSubset$expected <- c((rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[1], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1), 
                             (rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[2], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1))
    colnames(dataSubset) <- c("observed", "Class", "expected")
    lambdaGC <- round(c(GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[1]]), GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[2]])), digits=3)
    txt1 <- paste(unique(dataSubset$Class)[1], lambdaGC[1])
    txt2 <- paste(unique(dataSubset$Class)[2], lambdaGC[2])
    text <- paste("lambda gc", txt1, txt2, sep="\n")
    p1 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("Window-based: theta=", Thetas[thetaIndx])) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3.7, label=text, color="blue4", size=12)
    windowPlots[[thetaIndx]] <- p1
  }
  
  singlePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubset <- as.data.frame(ParameterEst[ParameterEst$theta==Thetas[thetaIndx], c(2,3)])
    dataSubset <- dataSubset[complete.cases(dataSubset), ]
    dataSubset$single <- as.numeric(dataSubset$single)
    dataSubset$expected <- c((rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[1], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1), 
                             (rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[2], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1))
    colnames(dataSubset) <- c("observed", "Class", "expected")
    lambdaGC <- round(c(GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[1]]), GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[2]])), digits=3)
    txt1 <- paste(unique(dataSubset$Class)[1], lambdaGC[1])
    txt2 <- paste(unique(dataSubset$Class)[2], lambdaGC[2])
    text <- paste("lambda gc", txt1, txt2, sep="\n")
    p2 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("theta=", Thetas[thetaIndx])) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3.7, label=text, color="blue4", size=15)
    singlePlots[[thetaIndx]] <- p2
  }
  
  combinePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubsetPre <- as.matrix(ParameterEst[ParameterEst$theta==Thetas[thetaIndx], c(1,2,3)])
    dataSubset <- rbind(dataSubsetPre[ , c(1,3)], dataSubsetPre[ , c(2,3)])
    dataSubset <- as.data.frame(dataSubset[complete.cases(dataSubset), ])
    colnames(dataSubset) <- c("observed", "Class")
    dataSubset$observed <- as.numeric(dataSubset$observed)
    dataSubset <- dataSubset %>% arrange(Class)
    dataSubset$expected <- c((rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[1], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1), 
                             (rank(dataSubset[dataSubset$Class==unique(dataSubset$Class)[2], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1))
    lambdaGC <- round(c(GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[1]]), GIF(dataSubset$observed[dataSubset$Class==unique(dataSubset$Class)[2]])), digits=3)
    txt1 <- paste(unique(dataSubset$Class)[1], lambdaGC[1])
    txt2 <- paste(unique(dataSubset$Class)[2], lambdaGC[2])
    text <- paste("lambda gc", txt1, txt2, sep="\n")
    p3 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("Combined: theta=", Thetas[thetaIndx])) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3.7, label=text, color="blue4", size=12)
    combinePlots[[thetaIndx]] <- p3
  }
  
  sixPanel <- do.call(ggarrange, c(c(windowPlots,singlePlots,combinePlots)[as.vector(sapply(1:length(Thetas), function(x) {c(x, length(Thetas)+x, 2*length(Thetas)+x)}))], nrow=length(Thetas), ncol=3, common.legend=TRUE, legend="bottom"))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related samples QQ plot: Variance=", VarianceValue, " & p=1000 & n=10000 & SKAT.haplotypes"), face="bold", size=40))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Tune h2: FDP inflation of related samples & P.causal=0.4 & rep=500 & MAC>=25 & n=5000"), face="bold", size=40))
  #png(paste0("/Users/xqi/Downloads/BinaryRelatedH0qqPlotsOld.png"), width=3500, height=3550)
  #png(paste0("/Users/xqi/Downloads/rareRelatedH0qqPlotsOldVar=",VarianceValue,".png"), width=3500, height=3550)
  #png(paste0("FDPInflatePlan2.png"), width=2500, height=2550)
  #png(paste0("FDPInflatePlan2.png"), width=2500, height=3200)
  plot(sixPanel)
  dev.off()
  
  return(list(H0Window=windowPlots,
              H0Single=singlePlots,
              H0Combine=combinePlots))
}

adjust <- parameterEstimateComparison(Thetas=c(0,4,7), 1)
H0Window <- adjust$H0Window
H0Combine <- adjust$H0Combine
H0Single <- adjust$H0Single



############# related H1 ################
# a function to rename methods
renameMethods <- function(originalName, TypeOfTest){
  if (originalName == "HMMK_ACAT"){
    newName <- "HMM knockoff, ACAT"
  } else if (originalName == "STAAR_identity"){
    newName <- "SCIT knockoffs, ACAT"
  } else if (originalName == "NormalK_ACAT"){
    newName <- 'SecondOrder knockoff, ACAT'
  } else if (originalName == "STAAR_identity_"){
    newName <- 'SCIP knockoff, ACAT'
  } else if (originalName == "STAAR"){
    if(TypeOfTest == "single variant") {
      newName <- "SCIT knockoffs, mixed model score test"
    } else {
      newName <- "SCIT knockoffs, STAAR" 
    }
  } else if (originalName == "HMMSK"){
    newName <- "HMM knockoff, score test"
  } else if (originalName == "MSK") {
    newName <- "SCIT knockoffs, score test"
  } else if (originalName == "SK") {
    newName <- "SCIP knockoff, score test"
  } else if (originalName == "NormalSK") {
    newName <- "SecondOrder knockoff, score test"
  } else if (originalName == "SummaryStat_score") {
    newName <- "SummaryStat-based knockoffs, score test"
  } else if (originalName == "SummaryStat_ACAT") {
    newName <- "SummaryStat-based knockoffs, ACAT"
  } else if (originalName == "SummaryStat_STAAR") {
    newName <- "SummaryStat-based knockoffs, mixed model score test"
  } else if (originalName == "Ghostknockoff_STAAR") {
    newName <- "Ghost knockoffs, mixed model score test"
  } else if (originalName ==  "Ghostknockoff_score") {
    newName <- "Ghost knockoffs, score test"
  } else if (originalName == "NormalK_multipleGLM") {
    newName <- "Multiple second order knockoffs, score test"
  } else if (originalName == "NormalK_multipleGLMM") {
    newName <- "Multiple second order knockoffs, mixed model score test"
  } else if (originalName == "NormalK_multipleACAT") {
    newName <- "Multiple second order knockoffs, ACAT"
  } else if (originalName == "NormalK_multipleSTAAR") {
    newName <- "Multiple second order knockoffs, STAAR"
  } else {
    newName <- originalName
  }
  return(newName)
}


combineResult_knockoff_rare_C502_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/combineResult_knockoff_C_23_10_0_detail.txt")
combineResult_knockoff_rare_C542_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/combineResult_knockoff_C_23_10_4_detail.txt")
combineResult_knockoff_rare_C582_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/combineResult_knockoff_C_23_10_7_detail.txt")

result_knockoff_rare_C502_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/result_knockoff_C_23_10_0_detail.txt")
result_knockoff_rare_C542_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/result_knockoff_C_23_10_4_detail.txt")
result_knockoff_rare_C582_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/fake/result_knockoff_C_23_10_7_detail.txt")

singleResult_knockoff_rare_C502_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=0/N=10/h2=1/singleResult_knockoff_C_1000_10_0_detail.txt")
singleResult_knockoff_rare_C542_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=4/N=10/h2=1/singleResult_knockoff_C_1000_10_4_detail.txt")
singleResult_knockoff_rare_C582_detail <- read.delim("/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/theta=7/N=10/h2=1/singleResult_knockoff_C_1000_10_7_detail.txt")

singleResult_knockoff_rare_C502_detail <- read.delim("/Users/xqi/Downloads/theta=0/singleResult_knockoff_D_737_10_0_detail.txt")
singleResult_knockoff_rare_C542_detail <- read.delim("/Users/xqi/Downloads/theta=4/singleResult_knockoff_D_737_10_4_detail.txt")
singleResult_knockoff_rare_C582_detail <- read.delim("/Users/xqi/Downloads/theta=7/singleResult_knockoff_D_737_10_7_detail.txt")

result02 <- list(singleResult.detail=singleResult_knockoff_rare_C502_detail,
                 result.detail=result_knockoff_rare_C502_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C502_detail)
result42 <- list(singleResult.detail=singleResult_knockoff_rare_C542_detail,
                 result.detail=result_knockoff_rare_C542_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C542_detail)
result82 <- list(singleResult.detail=singleResult_knockoff_rare_C582_detail,
                 result.detail=result_knockoff_rare_C582_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C582_detail)


#Thetas=c(0,4,7.5);EpsSDs=c(8,4,0.5);H0Result=H0Single;H2s=c(0.09,0.09,0.09,0.09)
#PCausal <- 0.1
#Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail)
#methodsForCompare=c(2,3,6,9,11,12,13,15,16);TypeOfTest="single variant"

# combine the above plots
# a function to combine 6-panel plot with 3-panel plot for window-based/single variant/combined knockoff-based inference methods comparison
compareRandomEffectsH1H0 <- function(Thetas, EpsSDs, Datasets, methodsForCompare, TypeOfTest, H0Result, H2s, PCausal) {
  # a matrix to store power estimates
  result.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(Thetas), ncol=8)
  for (thetaIndx in 1:length(Thetas)) {
    result.detail.power[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), c(7,8)] <- cbind(rep(Thetas[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)), rep(EpsSDs[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.power[(thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                              rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                              rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                                   rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                                   rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.power) <- c(colnames(Datasets)[1:4], "power", "Type", "theta", "epsSD")
  result.detail.power <- as.data.frame(result.detail.power)
  result.detail.power$Type <- as.factor(result.detail.power$Type)
  
  # a matrix to store fdp estimates
  result.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(Thetas), ncol=8)
  for (thetaIndx in 1:length(Thetas)) {
    result.detail.fdp[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), c(7,8)] <- cbind(rep(Thetas[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)), rep(EpsSDs[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.fdp[(thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                              rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                              rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                                   rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                                   rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.fdp) <- c(colnames(Datasets)[1:4], "fdp", "Type", "theta", "epsSD")
  result.detail.fdp <- as.data.frame(result.detail.fdp)
  result.detail.fdp$Type <- as.factor(result.detail.fdp$Type)
  
  # simplify Type labels
  result.detail.power$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.power$Type, "_median"), "single_"), "_combine"), "_power")
  result.detail.fdp$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.fdp$Type, "_median"), "single_"), "_combine"), "_fdp")
  # rename methods
  result.detail.power$Type <- unlist(sapply(result.detail.power$Type, function(x) renameMethods(x, TypeOfTest)))
  result.detail.fdp$Type <- unlist(sapply(result.detail.fdp$Type, function(x) renameMethods(x, TypeOfTest)))

  # combined data
  colnames(result.detail.power) <- colnames(result.detail.fdp) <- c(colnames(result.detail.fdp)[1:4], "value", colnames(result.detail.fdp)[6:8])
  combinedData <- rbind.data.frame(result.detail.fdp, result.detail.power)
  combinedData$valueType <- c(rep("FDP", times=nrow(result.detail.fdp)), rep("Power", times=nrow(result.detail.power)))
  
  # plot power comparison
  powerPlots <- fdpPlots <- combinePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    p1 <- ggplot(result.detail.power[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), ], aes(x=as.numeric(fdr), y=as.numeric(power), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 30), plot.title = element_text(size = 35)) +
      #ggtitle(paste0("Power comparison: theta=", Thetas[thetaIndx], " & Var(e)=", (EpsSDs[thetaIndx]), " h2= ", H2s[thetaIndx])) +
      ggtitle(paste0("Power comparison: theta=", Thetas[thetaIndx])) +
      xlab("False discovery rate") + ylab("Power") + 
      labs(fill = "Power type") +
      ylim(0, 1) 
    powerPlots[[thetaIndx]] <- p1
    
    p2 <- ggplot(result.detail.fdp[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), ], aes(x=as.numeric(fdr), y=as.numeric(fdp), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 30), plot.title = element_text(size = 35)) +
      ggtitle(paste0("FDP comparison: theta=", Thetas[thetaIndx])) +
      #ggtitle(paste0("FDP comparison: theta=", Thetas[thetaIndx], " & Var(e)=", (EpsSDs[thetaIndx]), " h2= ", H2s[thetaIndx])) +
      xlab("False discovery rate") + ylab("Observed FDR") + 
      labs(fill = "FDP type") +
      ylim(0, 0.4) +
      geom_abline(intercept=0, slope=1, linetype='dotted', col='black', size=1.5)
    fdpPlots[[thetaIndx]] <- p2
    
    thetaTemp <- ((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare))
    combinePlots[[thetaIndx]] <- ggplot(combinedData[c(thetaTemp, thetaTemp+nrow(result.detail.fdp)), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      facet_grid(valueType ~ ., scales="free_y") + 
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20), plot.title = element_text(size = 25)) +
      ggtitle(paste0("theta=", Thetas[thetaIndx])) +
      xlab("False discovery rate") + ylab("") +
      geom_abline(intercept=0, slope=1, linetype='dotted', col='black', size=2)
  }
  
  sixPanel <- do.call(ggarrange, c(combinePlots, nrow=3, ncol=1, common.legend=TRUE, legend="bottom"))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related samples: binary & ", TypeOfTest, " & 867 rep & p=1000 & h2=", unique(H2s), " & N.causal=", PCausal, " & SKAT.haplotypes & +- & cases from case fam & random control fam"), face="bold", size=40))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related binary outcomes & all samples"), face="bold", size=30))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related samples: binary & ", TypeOfTest, " & 867 rep & p=1000 & h2=", unique(H2s), " & N.causal=", PCausal, " & SKAT.haplotypes & +- & family as sample unit"), face="bold", size=40))
  # QQ plots under H0
  threePanel <- do.call(ggarrange, c(H0Result, nrow=3, ncol=1, common.legend=TRUE, legend="bottom"))
  threePanel <- annotate_figure(threePanel, top=text_grob(paste0("QQ plot of ", TypeOfTest, " under H0"), face="bold", size=30))
  
  overallPlot <- ggarrange(sixPanel, threePanel, nrow=1, ncol=2, widths=c(2, 1))
  #png(paste0("/Users/xqi/Downloads/Binary",TypeOfTest, "_", unique(H2s), "_", PCausal, "Overall.png"), width=4200, height=4200)
  png(paste0("/Users/xqi/Downloads/binary",TypeOfTest, "_", unique(H2s), "_", PCausal, "Overall.png"), width=2000, height=3000)
  plot(overallPlot)
  dev.off()
}

compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$result.detail,result42$result.detail,result82$result.detail),methodsForCompare=c(15,17,18,19,20),TypeOfTest="window-based",H0Result=H0Window,H2s=c(0.08,0.08,0.08,0.08),PCausal=0.4)
compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$combinedResult.detail,result42$combinedResult.detail,result82$combinedResult.detail),methodsForCompare=c(15,17,18,19,20),TypeOfTest="combined",H0Result=H0Combine,H2s=c(0.08,0.08,0.08,0.08),PCausal=0.4)
compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(2,3,6,9,11,12,13,15,16),TypeOfTest="single variant",H0Result=H0Single,H2s=c(0.08,0.08,0.08,0.08),PCausal=0.4)

compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$result.detail,result42$result.detail,result82$result.detail),methodsForCompare=c(15,17,22,24),TypeOfTest="window-based",H0Result=H0Window,H2s=0.5,PCausal=10)
compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$combinedResult.detail,result42$combinedResult.detail,result82$combinedResult.detail),methodsForCompare=c(15,17,22,24),TypeOfTest="combined",H0Result=H0Combine,H2s=1,PCausal=10)
compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),TypeOfTest="single variant",H0Result=H0Single,H2s=0.5,PCausal=10)

compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$result.detail,result42$result.detail,result82$result.detail),methodsForCompare=c(15,17,22,24),TypeOfTest="window-based",H0Result=H0Window,H2s=8,PCausal=10)
compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$combinedResult.detail,result42$combinedResult.detail,result82$combinedResult.detail),methodsForCompare=c(15,17,22,24),TypeOfTest="combined",H0Result=H0Combine,H2s=5,PCausal=10)



compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),TypeOfTest="single variant",H0Result=H0Single,H2s=2.5,PCausal=10)



compareRandomEffectsH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(2,11,13,15,16,17,19,20),TypeOfTest="single variant",H0Result=H0Single,H2s=5,PCausal=10)


###################### independent H0 #########################
ParameterEst <- aggregatePValue(NA, "/Users/xqi/Downloads/C/qqPlot/indep")
ParameterEst <- aggregatePValue(NA, "/Users/xqi/Downloads/qqPlot/indep")
colnames(ParameterEst) <- c("window", "single", "Class", "theta")

parameterEstimateComparison <- function(Thetas, VarianceValue) {
  windowPlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubset <- as.data.frame(ParameterEst[is.na(ParameterEst$theta), c(1,3)])
    dataSubset <- dataSubset[complete.cases(dataSubset), ]
    dataSubset$window <- as.numeric(dataSubset$window)
    dataSubset$expected <- (rank(dataSubset$window, ties.method="first")+.5)/(nrow(dataSubset)+1)
    colnames(dataSubset) <- c("observed", "Class", "expected")
    lambdaGC <- round(GIF(dataSubset$observed), digits=3)
    txt1 <- paste(unique(dataSubset$Class), lambdaGC)
    text <- paste("lambda gc", txt1, sep="\n")
    p1 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("Window-based: independent & theta=", Thetas[thetaIndx])) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3, label=text, color="blue4", size=12)
    windowPlots[[thetaIndx]] <- p1
  }
  
  singlePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubset <- as.data.frame(ParameterEst[is.na(ParameterEst$theta), c(2,3)])
    dataSubset <- dataSubset[complete.cases(dataSubset), ]
    dataSubset$single <- as.numeric(dataSubset$single)
    dataSubset$expected <- (rank(dataSubset$single, ties.method="first")+.5)/(nrow(dataSubset)+1)
    colnames(dataSubset) <- c("observed", "Class", "expected")
    lambdaGC <- round(GIF(dataSubset$observed), digits=3)
    txt1 <- paste(unique(dataSubset$Class), lambdaGC)
    text <- paste("lambda gc", txt1, sep="\n")
    p2 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("independent")) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3, label=text, color="blue4", size=15)
    singlePlots[[thetaIndx]] <- p2
  }
  
  combinePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubsetPre <- as.matrix(ParameterEst[is.na(ParameterEst$theta), c(1,2,3)])
    dataSubset <- rbind(dataSubsetPre[ , c(1,3)], dataSubsetPre[ , c(2,3)])
    dataSubset <- as.data.frame(dataSubset[complete.cases(dataSubset), ])
    colnames(dataSubset) <- c("observed", "Class")
    dataSubset$observed <- as.numeric(dataSubset$observed)
    dataSubset$expected <- (rank(dataSubset$observed, ties.method="first")+.5)/(nrow(dataSubset)+1)
    lambdaGC <- round(GIF(dataSubset$observed), digits=3)
    txt1 <- paste(unique(dataSubset$Class), lambdaGC)
    text <- paste("lambda gc", txt1, sep="\n")
    p3 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=6, aes(shape=Class, color=Class)) +
      theme_base() +
      scale_color_manual(values=c("#00AFBB", "#FC4E07")) +
      geom_abline(intercept=0, slope=1, col='#E69F00', size=2, linetype="twodash") +
      theme(legend.key=element_blank(), text=element_text(size=30), legend.text=element_text(size=30), legend.title=element_text(size=30), plot.title=element_text(size=30)) +
      ggtitle(paste0("Combined: independent & theta=", Thetas[thetaIndx])) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      ggplot2::annotate("text", 1, 3, label=text, color="blue4", size=12)
    combinePlots[[thetaIndx]] <- p3
  }
  
  sixPanel <- do.call(ggarrange, c(c(windowPlots,singlePlots,combinePlots)[as.vector(sapply(1:length(Thetas), function(x) {c(x, length(Thetas)+x, 2*length(Thetas)+x)}))], nrow=length(Thetas), ncol=3, common.legend=TRUE, legend="bottom"))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Independent samples QQ plot: binary & prevalence=0.1 & subsample CC ratio=1:1 & p=1000 & Variance=", VarianceValue), face="bold", size=40))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Independent samples QQ plot: keep mu-mean(mu) & p=1000 & n=4000"), face="bold", size=40))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Tune h2: FDP inflation of related samples & P.causal=0.4 & rep=500 & MAC>=25 & n=5000"), face="bold", size=40))
  #png(paste0("/Users/xqi/Downloads/BinaryIndependentH0qqPlots.png"), width=3500, height=1200)
  png(paste0("/Users/xqi/Downloads/BinaryIndependentH0qqPlotsVar=", VarianceValue,".png"), width=3500, height=1200)
  plot(sixPanel)
  dev.off()
  
  return(list(H0Window=windowPlots,
              H0Single=singlePlots,
              H0Combine=combinePlots))
}

adjust <- parameterEstimateComparison(Thetas=c(NA),1)
H0Window <- adjust$H0Window
H0Combine <- adjust$H0Combine
H0Single <- adjust$H0Single



###################### independent H1 #########################
# a function to rename methods
renameIndepMethods <- function(originalName, TypeOfTest){
  if (originalName == "HMMK_ACAT"){
    newName <- "HMM knockoff, ACAT"
  } else if (originalName == "STAAR_identity"){
    newName <- "SCIT knockoffs, ACAT"
  } else if (originalName == "NormalK_ACAT"){
    newName <- 'SecondOrder knockoff, ACAT'
  } else if (originalName == "STAAR_identity_"){
    newName <- 'SCIP knockoff, ACAT'
  } else if (originalName == "STAAR"){
    if(TypeOfTest == "single variant") {
      newName <- "SCIT knockoffs, mixed model score test"
    } else {
      newName <- "SCIT knockoffs, STAAR" 
    }
  } else if (originalName == "HMMSK"){
    newName <- "HMM knockoff, score test"
  } else if (originalName == "MSK") {
    newName <- "SCIT knockoffs, score test"
  } else if (originalName == "SK") {
    newName <- "SCIP knockoff, score test"
  } else if (originalName == "NormalSK") {
    newName <- "SecondOrder knockoff, score test"
  } else if (originalName == "SummaryStat_score") {
    newName <- "SummaryStat-based knockoffs, score test"
  } else if (originalName == "SummaryStat_ACAT") {
    newName <- "SummaryStat-based knockoffs, ACAT"
  } else if (originalName == "Ghostknockoff_STAAR") {
    newName <- "Ghost knockoffs, mixed model score test"
  } else if (originalName == "Ghostknockoff_ACAT") {
    newName <- "Ghost knockoffs, ACAT score test"
  } else if (originalName ==  "Ghostknockoff_score") {
    newName <- "Ghost knockoffs, score test"
  } else if (originalName ==  "Ghostknockoff_Wald") {
    newName <- "Ghost knockoffs, Wald test"
  } else if (originalName ==  "Ghostknockoff_LRT") {
    newName <- "Ghost knockoffs, likelihood ratio test"
  } else if (originalName == "NormalK_multipleGLM") {
    newName <- "Multiple second order knockoffs, score test"
  } else if (originalName == "NormalK_multipleACAT") {
    newName <- "Multiple second order knockoffs, ACAT"
  } else if (originalName == "NormalK_multipleGLMM") {
    newName <- "Multiple second order knockoffs, mixed model score test"
  } else {
    newName <- originalName
  }
  return(newName)
}
combineResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepCombineResult_knockoff_C_1000_10_NA_detail.txt")
singleResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepSingleResult_knockoff_C_1000_10_NA_detail.txt")
result_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepResult_knockoff_C_1000_10_NA_detail.txt")


combineResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/indepCombineResult_knockoff_D_1000_10_NA_detail.txt")
singleResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/IndepSingleResult_knockoff_D_1000_10_NA_detail.txt")
result_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/indepResult_knockoff_D_1000_10_NA_detail.txt")
indep0.8 <- list(result.detail=result_knockoff_rare_C5_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C5_detail,
                 singleResult.detail=singleResult_knockoff_rare_C5_detail)

#H2s=c(0.5);Datasets=as.data.frame(indep0.8$singleResult.detail);methodsForCompare=c(2,3,6,9,11,12,13);TypeOfTest="single variant";H0Result=H0Single
#H2s=c(0.08);PCausal=0.4;Datasets=as.data.frame(singleResult_knockoff_rare_C5_detail);methodsForCompare=c(2,3,6,9,11,12,13,14,15);TypeOfTest="single variant"
#H2s=c(0.08);PCausal=0.1;Datasets=as.data.frame(singleResult_knockoff_rare_C5_detail);methodsForCompare=c(2,4,5);TypeOfTest="single variant"
# a function to combine 6-panel plot with 3-panel plot for window-based/single variant/combined knockoff-based inference methods comparison
compareRandomEffectsH1H0 <- function(H2s, Datasets, methodsForCompare, TypeOfTest, H0Result, PCausal) {
  result.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.power[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.power[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(h2Indx-1)*nrow(Datasets)/length(H2s)), 
                              rep(TRUE, times=nrow(Datasets)/length(H2s)), 
                              rep(FALSE, times=(length(H2s)-h2Indx)*nrow(Datasets)/length(H2s))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(h2Indx-1)*nrow(Datasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(Datasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(Datasets)/length(H2s))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.power) <- c(colnames(Datasets)[1:4], "power", "Type", "hSquare")
  result.detail.power <- as.data.frame(result.detail.power)
  result.detail.power$fdr <- as.numeric(as.vector(result.detail.power$fdr))
  result.detail.power$power <- as.numeric(as.vector(result.detail.power$power))
  
  # a matrix to store fdp estimates
  result.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.fdp[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.fdp[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(h2Indx-1)*nrow(Datasets)/length(H2s)), 
                              rep(TRUE, times=nrow(Datasets)/length(H2s)), 
                              rep(FALSE, times=(length(H2s)-h2Indx)*nrow(Datasets)/length(H2s))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(h2Indx-1)*nrow(Datasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(Datasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(Datasets)/length(H2s))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.fdp) <- c(colnames(Datasets)[1:4], "fdp", "Type", "hSquare")
  result.detail.fdp <- as.data.frame(result.detail.fdp)
  result.detail.fdp$fdr <- as.numeric(as.vector(result.detail.fdp$fdr))
  result.detail.fdp$fdp <- as.numeric(as.vector(result.detail.fdp$fdp))
  
  # simplify Type labels
  result.detail.power$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.power$Type, "_median"), "single_"), "_combine"), "_power")
  result.detail.fdp$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.fdp$Type, "_median"), "single_"), "_combine"), "_fdp")
  # rename methods
  result.detail.power$Type <- unlist(sapply(result.detail.power$Type, function(x) renameIndepMethods(x, TypeOfTest)))
  result.detail.fdp$Type <- unlist(sapply(result.detail.fdp$Type, function(x) renameIndepMethods(x, TypeOfTest)))
  
  # combined data
  colnames(result.detail.power) <- colnames(result.detail.fdp) <- c(colnames(result.detail.fdp)[1:4], "value", colnames(result.detail.fdp)[6:7])
  combinedData <- rbind.data.frame(result.detail.fdp, result.detail.power)
  combinedData$valueType <- c(rep("FDP", times=nrow(result.detail.fdp)), rep("Power", times=nrow(result.detail.power)))
  

  # plot power comparison
  powerPlots <- fdpPlots <- combinePlots <- list()
  for (h2Indx in 1:length(H2s)) {
    p1 <- ggplot(result.detail.power[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), ], aes(x=as.numeric(fdr), y=as.numeric(power), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.03)) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 30), plot.title = element_text(size = 35)) +
      ggtitle(paste0("Power comparison")) +
      xlab("False discovery rate") + ylab("Power") + 
      labs(fill = "Power type") +
      ylim(0, 1) +
      theme(plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)) # Left margin
    powerPlots[[h2Indx]] <- p1
    
    p2 <- ggplot(result.detail.fdp[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), ], aes(x=as.numeric(fdr), y=as.numeric(fdp), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.03)) +
      # scale_x_discrete("False discovery rate", labels=c(sapply(as.character(seq(from=0, to=0.2, by=0.02)), function(x) c(x, "")))[1:length(seq(from=0, to=0.2, by=0.01))]) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 30), legend.text = element_text(size = 30), legend.title = element_text(size = 30), plot.title = element_text(size = 35)) +
      ggtitle(paste0("FDP comparison")) +
      xlab("False discovery rate") + ylab("Observed FDR") + 
      labs(fill = "FDP type") +
      ylim(0, 0.3) +
      geom_abline(intercept=0, slope=1, linetype='dotted', col='black', size=1.5) +
      theme(plot.margin = margin(t = 0,  # Top margin
                                 r = 0,  # Right margin
                                 b = 0,  # Bottom margin
                                 l = 0)) # Left margin
    fdpPlots[[h2Indx]] <- p2
    
    thetaTemp <- ((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(methodsForCompare))
    combinePlots[[h2Indx]] <- ggplot(combinedData[c(thetaTemp, thetaTemp+nrow(result.detail.fdp)), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
      geom_line(aes(color=Type), size=1.5) +
      geom_point(aes(color=Type), size=4) +
      facet_grid(valueType ~ ., scales="free_y") + 
      theme_base() +
      scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
      theme(legend.key = element_blank(), text = element_text(size = 20), legend.text = element_text(size = 20), legend.title = element_text(size = 20), plot.title = element_text(size = 25)) +
      ggtitle("independent") +
      xlab("False discovery rate") + ylab("") +
      geom_abline(intercept=0, slope=1, linetype='dotted', col='black', size=2)
  }
  
  sixPanel <- do.call(ggarrange, c(combinePlots, nrow=1, ncol=1, common.legend=TRUE, legend="bottom"))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Independent binary & h2=", H2s, " & N.causal=", PCausal), face="bold", size=30))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Independent samples: continuous & ", TypeOfTest, " & 347 rep & p=1000 & h2=", H2s, " & N.causal=", PCausal, " & SKAT.haplotypes & +- & indiv as sample unit"), face="bold", size=40))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Independent samples: ", TypeOfTest, " & 100 rep & n=4000 & p=1000 & h2=", H2s, " & P.causal=", PCausal, " & 3 phased files & empirical cor"), face="bold", size=40))
  #png(paste0("/Users/xqi/Downloads/empiricalCorIndep", TypeOfTest, "_", H2s, "_", PCausal, "Overall.png"), width=2500, height=1300)
  #plot(sixPanel)
  #dev.off()
  threePanel <- do.call(ggarrange, c(H0Result, nrow=1, ncol=1, common.legend=TRUE, legend="bottom"))
  threePanel <- annotate_figure(threePanel, top=text_grob(paste0("QQ plot of ", TypeOfTest, " under H0"), face="bold", size=30))
  
  overallPlot <- ggarrange(sixPanel, threePanel, nrow=1, ncol=2, widths=c(2, 1))
  #png(paste0("/Users/xqi/Downloads/BinaryIndep", TypeOfTest, "_", H2s, "_", PCausal, "Overall.png"), width=4200, height=1430)
  png(paste0("/Users/xqi/Downloads/binaryIndep", TypeOfTest, "_", H2s, "_", PCausal, "Overall.png"), width=2000, height=1000)
  plot(overallPlot)
  dev.off()
}

compareRandomEffectsH1H0(H2s=c(0.08),PCausal=0.4,Datasets=as.data.frame(indep0.8$result.detail),methodsForCompare=c(15,16,17,18),TypeOfTest="window-based",H0Result=H0Window)
compareRandomEffectsH1H0(H2s=c(0.08),PCausal=0.4,Datasets=as.data.frame(indep0.8$combinedResult.detail),methodsForCompare=c(15,16,17,18),TypeOfTest="combined",H0Result=H0Combine)
compareRandomEffectsH1H0(H2s=c(0.08),PCausal=0.4,Datasets=as.data.frame(indep0.8$singleResult.detail),methodsForCompare=c(2,3,6,9,11,12,13,14,15),TypeOfTest="single variant",H0Result=H0Single)

compareRandomEffectsH1H0(H2s=c(0.5),PCausal=10,Datasets=as.data.frame(indep0.8$result.detail),methodsForCompare=c(15,20),TypeOfTest="window-based",H0Result=H0Window)
compareRandomEffectsH1H0(H2s=c(1),PCausal=10,Datasets=as.data.frame(indep0.8$combinedResult.detail),methodsForCompare=c(15,20),TypeOfTest="combined",H0Result=H0Combine)
compareRandomEffectsH1H0(H2s=c(2.5),PCausal=10,Datasets=as.data.frame(indep0.8$singleResult.detail),methodsForCompare=c(13,16),TypeOfTest="single variant",H0Result=H0Single)

compareRandomEffectsH1H0(H2s=c(8),PCausal=10,Datasets=as.data.frame(indep0.8$result.detail),methodsForCompare=c(15,20),TypeOfTest="window-based",H0Result=H0Window)
compareRandomEffectsH1H0(H2s=c(5),PCausal=10,Datasets=as.data.frame(indep0.8$combinedResult.detail),methodsForCompare=c(15,20),TypeOfTest="combined",H0Result=H0Combine)



compareRandomEffectsH1H0(H2s=c(2.5),PCausal=10,Datasets=as.data.frame(indep0.8$singleResult.detail),methodsForCompare=c(13,16,18,19),TypeOfTest="single variant",H0Result=H0Single)




combineResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/indepCombineResult_knockoff_D_100_10_NA_detail.txt")
singleResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/indepSingleResult_knockoff_D_100_10_NA_detail.txt")
result_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/indep/indepResult_knockoff_D_100_10_NA_detail.txt")
indep0.8 <- list(result.detail=result_knockoff_rare_C5_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C5_detail,
                 singleResult.detail=singleResult_knockoff_rare_C5_detail)
compareRandomEffectsH1H0(H2s=c(2.77),PCausal=10,Datasets=as.data.frame(indep0.8$result.detail),methodsForCompare=c(15,20),TypeOfTest="window-based",H0Result=H0Window)
compareRandomEffectsH1H0(H2s=c(2.77),PCausal=10,Datasets=as.data.frame(indep0.8$singleResult.detail),methodsForCompare=c(13,16),TypeOfTest="single variant",H0Result=H0Single)



combineResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepCombineResult_knockoff_C_100_10_NA_detail.txt")
singleResult_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepSingleResult_knockoff_C_100_10_NA_detail.txt")
result_knockoff_rare_C5_detail <- read.delim("/Users/xqi/Downloads/C/indep/indepResult_knockoff_C_100_10_NA_detail.txt")
indep0.8 <- list(result.detail=result_knockoff_rare_C5_detail,
                 combinedResult.detail=combineResult_knockoff_rare_C5_detail,
                 singleResult.detail=singleResult_knockoff_rare_C5_detail)
compareRandomEffectsH1H0(H2s=c(1),PCausal=10,Datasets=as.data.frame(indep0.8$result.detail),methodsForCompare=c(15,20),TypeOfTest="window-based",H0Result=H0Window)
compareRandomEffectsH1H0(H2s=c(1),PCausal=10,Datasets=as.data.frame(indep0.8$singleResult.detail),methodsForCompare=c(2,13,16),TypeOfTest="single variant",H0Result=H0Single)




barPlots <- function(Thetas, EpsSDs, Datasets, methodsForCompare, indepDatasets, indepmethodsForCompare, TypeOfTest, H2s, PCausal, FDPSpecify) {
  # a matrix to store power estimates
  result.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(Thetas), ncol=8)
  for (thetaIndx in 1:length(Thetas)) {
    result.detail.power[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), c(7,8)] <- cbind(rep(Thetas[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)), rep(EpsSDs[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.power[(thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                              rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                              rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                                   rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                                   rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.power) <- c(colnames(Datasets)[1:4], "power", "Type", "theta", "epsSD")
  result.detail.power <- as.data.frame(result.detail.power)
  result.detail.power$Type <- as.factor(result.detail.power$Type)
  
  # a matrix to store fdp estimates
  result.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(methodsForCompare)*length(Thetas), ncol=8)
  for (thetaIndx in 1:length(Thetas)) {
    result.detail.fdp[((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare)), c(7,8)] <- cbind(rep(Thetas[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)), rep(EpsSDs[thetaIndx], times=length(seq(0,0.2,by=0.01))*length(methodsForCompare)))
    for (replicateIndx in 1:length(methodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.fdp[(thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                              rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                              rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                              Datasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(Datasets[c(rep(FALSE, times=(thetaIndx-1)*nrow(Datasets)/length(Thetas)), 
                                   rep(TRUE, times=nrow(Datasets)/length(Thetas)), 
                                   rep(FALSE, times=(length(Thetas)-thetaIndx)*nrow(Datasets)/length(Thetas))) & 
                                   Datasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(Datasets)[4+(ncol(Datasets)-4)/2+methodsForCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.fdp) <- c(colnames(Datasets)[1:4], "fdp", "Type", "theta", "epsSD")
  result.detail.fdp <- as.data.frame(result.detail.fdp)
  result.detail.fdp$Type <- as.factor(result.detail.fdp$Type)
  
  # simplify Type labels
  result.detail.power$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.power$Type, "_median"), "single_"), "_combine"), "_power")
  result.detail.fdp$Type <- str_remove(str_remove(str_remove(str_remove(result.detail.fdp$Type, "_median"), "single_"), "_combine"), "_fdp")
  # rename methods
  result.detail.power$Type <- unlist(sapply(result.detail.power$Type, function(x) renameMethods(x, TypeOfTest)))
  result.detail.fdp$Type <- unlist(sapply(result.detail.fdp$Type, function(x) renameMethods(x, TypeOfTest)))
  
  indepresult.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    indepresult.detail.power[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare))
    for (replicateIndx in 1:length(indepmethodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        indepresult.detail.power[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                   indepDatasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                        rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                        rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                        indepDatasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+indepmethodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasets)[4+indepmethodsForCompare[replicateIndx]]))
      }
    }
  }
  indepresult.detail.power <- cbind.data.frame(indepresult.detail.power, "indep")
  colnames(indepresult.detail.power) <- c(colnames(indepDatasets)[1:4], "power", "Type", "hSquare", "theta")
  indepresult.detail.power <- rbind.data.frame(indepresult.detail.power, indepresult.detail.power)
  indepresult.detail.power$Type <- c(indepresult.detail.power$Type[1:(nrow(indepresult.detail.power)/2)],
                                     rep(ifelse(unique(indepresult.detail.power$Type)[1]=="NormalK_median_multipleGLM_power", "NormalK_median_multipleGLMM_power", "Ghostknockoff_STAAR_power"), times=nrow(indepresult.detail.power)/4),
                                     rep(ifelse(unique(indepresult.detail.power$Type)[2]=="Ghostknockoff_score_power", "Ghostknockoff_STAAR_power", "NormalK_median_multipleGLMM_power"), times=nrow(indepresult.detail.power)/4))
  
  # a matrix to store fdp estimates
  indepresult.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    indepresult.detail.fdp[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare))
    for (replicateIndx in 1:length(indepmethodsForCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        indepresult.detail.fdp[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepmethodsForCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                   indepDatasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                        rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                        rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                        indepDatasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(indepDatasets)-4)/2+indepmethodsForCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasets)[4+(ncol(indepDatasets)-4)/2+indepmethodsForCompare[replicateIndx]]))
      }
    }
  }
  indepresult.detail.fdp <- cbind.data.frame(indepresult.detail.fdp, "indep")
  colnames(indepresult.detail.fdp) <- c(colnames(indepDatasets)[1:4], "fdp", "Type", "hSquare", "theta")
  indepresult.detail.fdp <- rbind.data.frame(indepresult.detail.fdp, indepresult.detail.fdp)
  indepresult.detail.fdp$Type <- c(indepresult.detail.fdp$Type[1:(nrow(indepresult.detail.fdp)/2)],
                                   rep(ifelse(unique(indepresult.detail.fdp$Type)[1]=="NormalK_median_multipleGLM_fdp", "NormalK_median_multipleGLMM_fdp", "Ghostknockoff_STAAR_fdp"), times=nrow(indepresult.detail.fdp)/4),
                                   rep(ifelse(unique(indepresult.detail.fdp$Type)[2]=="Ghostknockoff_score_fdp", "Ghostknockoff_STAAR_fdp", "NormalK_median_multipleGLMM_fdp"), times=nrow(indepresult.detail.fdp)/4))
  
  # simplify Type labels
  indepresult.detail.power$Type <- str_remove(str_remove(str_remove(str_remove(indepresult.detail.power$Type, "_median"), "single_"), "_combine"), "_power")
  indepresult.detail.fdp$Type <- str_remove(str_remove(str_remove(str_remove(indepresult.detail.fdp$Type, "_median"), "single_"), "_combine"), "_fdp")
  # rename methods
  indepresult.detail.power$Type <- unlist(sapply(indepresult.detail.power$Type, function(x) renameIndepMethods(x, TypeOfTest)))
  indepresult.detail.fdp$Type <- unlist(sapply(indepresult.detail.fdp$Type, function(x) renameIndepMethods(x, TypeOfTest)))
  
  result.detail.power <- rbind.data.frame(indepresult.detail.power[,-(ncol(indepresult.detail.power)-1)], result.detail.power[,-ncol(result.detail.power)])
  result.detail.power <- result.detail.power[result.detail.power$fdr==FDPSpecify,]
  result.detail.power$power <- as.numeric(result.detail.power$power)
  result.detail.fdp <- rbind.data.frame(indepresult.detail.fdp[,-(ncol(indepresult.detail.fdp)-1)], result.detail.fdp[,-ncol(result.detail.fdp)])
  result.detail.fdp <- result.detail.fdp[result.detail.fdp$fdr==FDPSpecify,]
  result.detail.fdp$fdp <- as.numeric(result.detail.fdp$fdp)
  
  powerPlots <- ggplot(result.detail.power, aes(fill=Type, y=power, x=theta)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(name="Sample relatedness", labels=c("indep"="Independent", "0"="Related\ngenotypes", "4"="Related\ngeno & pheno\ntheta=4", "7"="Related\ngeno & pheno\ntheta=7"),limits=c("indep","0","4","7")) +
    scale_y_continuous(name="Observed power", breaks=seq(from=0, to=1, by=0.1)) +
    theme_base() +
    ggtitle(paste0("Power comparison")) +
    ylim(0, 1) +
    labs(fill = "Methods")
  
  fdpPlots <- ggplot(result.detail.fdp, aes(fill=Type, y=fdp, x=theta)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(name="Sample relatedness", labels=c("indep"="Independent", "0"="Related\ngenotypes", "4"="Related\ngeno & pheno\ntheta=4", "7"="Related\ngeno & pheno\ntheta=7"),limits=c("indep","0","4","7")) +
    scale_y_continuous(name="Observed FDR", breaks=seq(from=0, to=0.5, by=0.05)) +
    theme_base() +
    ggtitle(paste0("FDP comparison")) +
    labs(fill = "Methods") +
    geom_abline(intercept=FDPSpecify, slope=0, linetype='dotted', col='black', size=1.5)
  
  Plots <- list()
  Plots[[1]] <- powerPlots; Plots[[2]] <- fdpPlots 
  sixPanel <- do.call(ggarrange, c(Plots, nrow=1, ncol=2, common.legend=TRUE, legend="bottom"))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Binary outcome: FDR=", FDPSpecify, " & ", TypeOfTest, " & n=10000 & p=1000 & h2=", unique(H2s), " & N.causal=", PCausal), face="bold", size=22))
  #png(paste0("/Users/xqi/Downloads/Continuous",TypeOfTest, "_", unique(H2s), "_", PCausal, "Overall.png"), width=1200, height=660)
  #plot(sixPanel)
  #dev.off()
}

barPlots(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$result.detail,result42$result.detail,result82$result.detail),methodsForCompare=c(15,17,22,24),indepDatasets=as.data.frame(indep0.8$result.detail),indepmethodsForCompare=c(15,20),TypeOfTest="window-based",H2s=8,PCausal=10,FDPSpecify=0.1)
barPlots(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepmethodsForCompare=c(13,16),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.1)


Zdistribution <- read.delim("~/Downloads/relatedDensityZ_D_0.9_detail.txt")

ZPlots <- list()
for (colIndx in 1:ncol(Zdistribution)) {
  df <- as.data.frame(Zdistribution[,colIndx])
  colnames(df) <- "zValue"
  p <- ggplot(df, aes(x=zValue)) + 
    geom_density() +
    geom_vline(aes(xintercept=mean(zValue)), color="blue", linetype="dashed", size=1) +
    labs(title=paste0("theta=",c(0,4,7)[colIndx]) ,x="sigmoid(mu+b0_pheno)", y="Density")
  ZPlots[[colIndx]] <- p
}
sixPanel <- do.call(ggarrange, c(ZPlots, nrow=1, ncol=3, common.legend=TRUE, legend="bottom"))
sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Binary outcome: related samples & prevalence=0.9"), face="bold", size=22))
png(paste0("/Users/xqi/Downloads/ContinuousZScore_",0.9,"_Overall.png"), width=1200, height=480)
plot(sixPanel)
dev.off()
