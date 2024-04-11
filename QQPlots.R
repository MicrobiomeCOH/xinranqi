###################################################
### used for drawing QQ plots of four scenarios 
###################################################

library(dplyr)
library(ggplot2)
library("ggpubr")
library(patchwork)
library(stringr)
library(ggrepel)
library(gridExtra)
library(ggthemes)
library(grid)

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

ACATSTAAR_0.4_0_2_detail <- aggregatePValue(0, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeD/D/theta=0")
ACATSTAAR_0.4_4_2_detail <- aggregatePValue(4, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeD/D/theta=4")
ACATSTAAR_0.4_7_2_detail <- aggregatePValue(7, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeD/D/theta=7")
ParameterEst <- rbind.data.frame(ACATSTAAR_0.4_0_2_detail,  ACATSTAAR_0.4_4_2_detail, ACATSTAAR_0.4_7_2_detail)
colnames(ParameterEst) <- c("window", "single", "Methods", "theta")
ParameterEst <- ParameterEst %>% arrange(Methods)

###################### independent H0 #########################
indepParameterEst  <- aggregatePValue(NA, "/oak/stanford/groups/zihuai/XinranQi/Project1/continuous/C/indep")
indepParameterEst <- aggregatePValue(NA, "/oak/stanford/groups/zihuai/XinranQi/Project1/SchemeA/D/indep")
colnames(indepParameterEst) <- c("window", "single", "Methods", "theta")


parameterEstimateComparison <- function(Thetas) {
  # independent QQ plot
  dataSubset <- as.data.frame(indepParameterEst[is.na(indepParameterEst$theta), c(2,3)])
  dataSubset <- dataSubset[complete.cases(dataSubset), ]
  dataSubset$single <- as.numeric(dataSubset$single)
  dataSubset$expected <- (rank(dataSubset$single, ties.method="first")+.5)/(nrow(dataSubset)+1)
  colnames(dataSubset) <- c("observed", "Methods", "expected")
  dataSubset[dataSubset=="ACAT"] <- "Score test"
  dataSubset[dataSubset=="STAAR"] <- "Mixed model score test"
  
  lambdaGC <- round(GIF(dataSubset$observed), digits=3)
  txt1 <- paste(unique(dataSubset$Methods), lambdaGC)
  text <- paste("Lambda gc", txt1, sep="\n")
  annotations <- data.frame(xpos=-Inf,ypos=Inf,annotateText=text,hjustvar=-0.25,vjustvar=1.5) #<- adjust

  indepQQ <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
    geom_point(size=3, shape=17, aes(color=Methods)) +
    scale_color_manual(values=c("#016c59","#b30000")) +
    geom_abline(intercept=0, slope=1, col='grey', size=1.5) +
    #ggtitle("Independent genotypes & quantitative phenotypes") +
    ggtitle("Independent genotypes & dichotomous phenotypes") +
    xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
    geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=4.5) +
    #ggplot2::annotate("text", 1.5, 5, label=text, color="blue4") +
    theme_bw() +    
    theme(legend.position="bottom",
          plot.title=element_text(face="bold", size=12),
          axis.ticks=element_line(colour="grey70", size=0.2),
          panel.grid.major=element_line(colour="grey70", size=0.2),
          panel.grid.minor=element_blank())

  # related QQ plots
  singlePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    dataSubset <- as.data.frame(ParameterEst[ParameterEst$theta==Thetas[thetaIndx], c(2,3)])
    dataSubset <- dataSubset[complete.cases(dataSubset), ]
    dataSubset$single <- as.numeric(dataSubset$single)
    dataSubset$expected <- c((rank(dataSubset[dataSubset$Methods==unique(dataSubset$Methods)[1], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1), 
                             (rank(dataSubset[dataSubset$Methods==unique(dataSubset$Methods)[2], 1], ties.method="first")+.5)/(nrow(dataSubset)/2+1))
    colnames(dataSubset) <- c("observed", "Methods", "expected")
    dataSubset[dataSubset=="ACAT"] <- "Score test"
    dataSubset[dataSubset=="STAAR"] <- "Mixed model score test"
    
    lambdaGC <- round(c(GIF(dataSubset$observed[dataSubset$Methods==unique(dataSubset$Methods)[1]]), GIF(dataSubset$observed[dataSubset$Methods==unique(dataSubset$Methods)[2]])), digits=3)
    txt1 <- paste(unique(dataSubset$Methods)[1], lambdaGC[1])
    txt2 <- paste(unique(dataSubset$Methods)[2], lambdaGC[2])
    text <- paste("Lambda gc", txt1, txt2, sep="\n")
    annotations1 <- data.frame(xpos=-Inf,ypos=Inf,annotateText=text,hjustvar=-0.25,vjustvar=1.5) #<- adjust
    
    p2 <- ggplot(dataSubset, aes(x=-log10(expected), y=-log10(observed))) + 
      geom_point(size=3, aes(shape=Methods, color=Methods)) +
      scale_color_manual(values=c("#b30000", "#016c59")) +
      geom_abline(intercept=0, slope=1, col='grey', size=1.5) +
      ggtitle(ifelse(thetaIndx==1, "Related genotypes", paste0("Related geno & pheno \u03B8=", Thetas[thetaIndx]))) +
      xlab("Expected -log(p-value)") + ylab("Observed -log(p-value)") +
      geom_text(data=annotations1,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), size=4.5) +
      #ggplot2::annotate("text", 1.5, 5, label=text, color="blue4")  +
      labs(color="Methods") +
      theme_bw() +    
      theme(plot.title=element_text(face="bold", size=12),
            axis.ticks=element_line(colour="grey70", size=0.2),
            panel.grid.major=element_line(colour="grey70", size=0.2),
            panel.grid.minor=element_blank())
    singlePlots[[thetaIndx]] <- p2
  }

  sixPanel <- do.call(ggarrange, c(singlePlots, nrow=length(Thetas), ncol=1, common.legend=TRUE, legend="bottom"))
  #sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related quantitative phenotypes"), face="bold"))
  sixPanel <- annotate_figure(sixPanel, top= text_grob(paste0("Related dichotomous phenotypes; sampling scheme C"), face="bold"))
  sixPanel2 <- ggarrange(indepQQ, sixPanel, nrow=2, ncol=1, heights=c(1.05, 3))
  #c(1.05, 3)
  #c(1,3)
  
  return(list(H0Single=sixPanel2))
}

continuous <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single
PlanB <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single
PlanA <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single
#PlanD <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single
PlanC <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single

Plan11 <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single
Plan12 <- parameterEstimateComparison(Thetas=c(0,4,7))$H0Single


sixPanel3 <- ggarrange(continuous, Plan6, nrow=1, ncol=2)
sixPanel4 <- ggarrange(Plan9, Plan8, Plan7, nrow=1, ncol=3)

sixPanel5 <- ggarrange(Plan11, Plan12, nrow=1, ncol=2)
sixPanel3 <- ggarrange(continuous, PlanA, PlanB, PlanC, nrow=1, ncol=4)

sixPanel3 <- PlanC

png(paste0("/oak/stanford/groups/zihuai/XinranQi/Project1/QQFourScenarios1.png"), width=350*5, height=1500*1.25)
png(paste0("/oak/stanford/groups/zihuai/XinranQi/Project1/QQFourScenariosNew.png"), width=350*4, height=1500)
plot(sixPanel3)  
#plot(sixPanel4)
#plot(sixPanel5)  
dev.off()
