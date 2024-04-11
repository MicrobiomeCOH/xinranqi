compareH1H0 <- function(Thetas, EpsSDs, Datasets, methodsForCompare, indepDatasets, indepDatasetsMeta,
                        indepmethodsForCompare, indepWaldLRTCompare, indepWaldLRTCompareMeta, TypeOfTest, H2s, PCausal, FDPSpecify, Kest) {
  ### Part 1 Related samples
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
  
  # aggregated FDP & power curves
  combinePlots <- list()
  for (thetaIndx in 1:length(Thetas)) {
    thetaTemp <- ((thetaIndx-1)*length(seq(0,0.2,by=0.01))*length(methodsForCompare)+1):(thetaIndx*length(seq(0,0.2,by=0.01))*length(methodsForCompare))
    combinePlots[[thetaIndx]] <- ggplot(combinedData[c(thetaTemp, thetaTemp+nrow(result.detail.fdp)), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
      scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
      geom_line(aes(color=Type), size=2) +
      facet_grid_sc(rows=vars(valueType), scales=list(y=scales_y), switch="y", labeller=as_labeller(c(FDP="Observed FDR", Power="Power"))) + 
      geom_abline(data=data.frame(yint=0,yslope=1,valueType="FDP"), aes(intercept=yint, slope=yslope), col='grey', size=2) +
      ggtitle(ifelse(thetaIndx==1, "Related genotypes", paste0("Related geno & pheno \u03B8=", Thetas[thetaIndx]))) +
      scale_colour_brewer(type="seq", palette="Spectral") +
      labs(color="Methods") +
      ylab(NULL) +
      theme_bw() +    
      theme(plot.title = element_text(face = "bold", size = 12),
            axis.ticks = element_line(colour = "grey70", size = 0.2),
            panel.grid.major = element_line(colour = "grey70", size = 0.2),
            panel.grid.minor = element_blank())
  }


  ### Part 2 bar plot for related samples
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
  colnames(indepresult.detail.power) <- c(colnames(indepDatasets)[1:4], "value", "Type", "hSquare", "theta")
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
  colnames(indepresult.detail.fdp) <- c(colnames(indepDatasets)[1:4], "value", "Type", "hSquare", "theta")
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
  result.detail.power$value <- as.numeric(result.detail.power$value)
  result.detail.fdp <- rbind.data.frame(indepresult.detail.fdp[,-(ncol(indepresult.detail.fdp)-1)], result.detail.fdp[,-ncol(result.detail.fdp)])
  result.detail.fdp <- result.detail.fdp[result.detail.fdp$fdr==FDPSpecify,]
  result.detail.fdp$value <- as.numeric(result.detail.fdp$value)
  
  combinedData01 <- rbind.data.frame(result.detail.fdp, result.detail.power)
  combinedData01$valueType <- c(rep("FDP", times=nrow(result.detail.fdp)), rep("Power", times=nrow(result.detail.power)))
  
  # aggregated bar plots
  barCombinePlots <- ggplot(combinedData01, aes(fill=Type, y=value, x=theta)) + 
    geom_bar(position="dodge", stat="identity") +
    scale_x_discrete(name="Sample relatedness", labels=c("indep"="Independent", "0"="Related\ngenotypes", "4"="Related\ngeno & pheno\n\u03B8=4", "7"="Related\ngeno & pheno\n\u03B8=7"),limits=c("indep","0","4","7")) +
    ylab(NULL) +
    #ggtitle(paste0("Quantitative phenotypes")) +
    ggtitle(paste0("Dichotomous phenotypes")) +
    #ggtitle(paste0("Dichotomous phenotypes; sampling scheme C; K=", Kest)) +
    labs(fill="Methods") +
    facet_grid_sc(rows=vars(valueType), scales=list(y=scales_y), switch="y", labeller=as_labeller(c(FDP="Observed FDR", Power="Power"))) + 
    geom_hline(data=data.frame(yint=FDPSpecify,valueType="FDP"), aes(yintercept=yint), col='grey', size=2) +
    scale_fill_brewer(type="seq", palette="Spectral") +
    theme_bw() +    
    theme(plot.title = element_text(face = "bold", size = 12),
      #legend.background = element_rect(fill = "white", size = 4, colour = "white"),
      #legend.justification = c(0, 1),
      #legend.position = c(0, 1),
      axis.ticks = element_line(colour = "grey70", size = 0.2),
      panel.grid.major = element_line(colour = "grey70", size = 0.2),
      panel.grid.minor = element_blank()
    )

  
  ### Part 3 Independent samples
  result.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.power[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare))
    for (replicateIndx in 1:length(indepWaldLRTCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.power[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                   indepDatasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                        rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                        rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                        indepDatasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+indepWaldLRTCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasets)[4+indepWaldLRTCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.power) <- c(colnames(indepDatasets)[1:4], "power", "Type", "hSquare")
  result.detail.power <- as.data.frame(result.detail.power)
  result.detail.power$fdr <- as.numeric(as.vector(result.detail.power$fdr))
  result.detail.power$power <- as.numeric(as.vector(result.detail.power$power))
  
  # a matrix to store fdp estimates
  result.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.fdp[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare))
    for (replicateIndx in 1:length(indepWaldLRTCompare)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.fdp[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompare) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                   rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                   rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                   indepDatasets[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasets[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasets)/length(H2s)), 
                                        rep(TRUE, times=nrow(indepDatasets)/length(H2s)), 
                                        rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasets)/length(H2s))) & 
                                        indepDatasets[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(indepDatasets)-4)/2+indepWaldLRTCompare[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasets)[4+(ncol(indepDatasets)-4)/2+indepWaldLRTCompare[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.fdp) <- c(colnames(indepDatasets)[1:4], "fdp", "Type", "hSquare")
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
  
  indepCombinePlots <- ggplot(combinedData[combinedData$Type %in% c("GhostKnockoff, score test","GhostKnockoff, Wald test","GhostKnockoff, likelihood ratio test"), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
    scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
    geom_line(aes(color=Type), size=2) +
    facet_grid_sc(rows=vars(valueType), scales=list(y=scales_y), switch="y", labeller=as_labeller(c(FDP="Observed FDR", Power="Power"))) + 
    geom_abline(data=data.frame(yint=0,yslope=1,valueType="FDP"), aes(intercept=yint, slope=yslope), col='grey', size=2) +
    #ggtitle("Quantitative phenotypes") +
    ggtitle("Dichotomous phenotypes") +
    scale_color_manual(values=c("#6a51a3", "#fe9929", "#238b45")) +
    labs(color="Methods") +
    ylab(NULL) +
    theme_bw() +    
    theme(legend.position = c(0.37, 0.87),
          plot.title = element_text(face = "bold", size = 12),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  indepPlots <- ggplot(combinedData[combinedData$Type %in% c("GhostKnockoff, score test","IndividualData knockoff, score test"), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
    scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
    geom_line(aes(color=Type), size=2) +
    facet_grid_sc(rows=vars(valueType), scales=list(y=scales_y), switch="y", labeller=as_labeller(c(FDP="Observed FDR", Power="Power"))) + 
    geom_abline(data=data.frame(yint=0,yslope=1,valueType="FDP"), aes(intercept=yint, slope=yslope), col='grey', size=2) +
    #ggtitle("Independent genotypes & quantitative phenotypes") +
    ggtitle("Independent genotypes & dichotomous phenotypes") +
    scale_color_manual(values=c("#fe9929", "#3182bd")) +
    labs(color="Methods") +
    ylab(NULL) +
    theme_bw() +    
    theme(legend.position = c(0.4, 0.88),
          plot.title = element_text(face = "bold", size = 12),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())
  
  
  ### Part 4 Independent samples: meta analysis
  result.detail.power <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.power[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta))
    for (replicateIndx in 1:length(indepWaldLRTCompareMeta)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.power[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasetsMeta[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasetsMeta)/length(H2s)), 
                                       rep(TRUE, times=nrow(indepDatasetsMeta)/length(H2s)), 
                                       rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasetsMeta)/length(H2s))) & 
                                       indepDatasetsMeta[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasetsMeta[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasetsMeta)/length(H2s)), 
                                        rep(TRUE, times=nrow(indepDatasetsMeta)/length(H2s)), 
                                        rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasetsMeta)/length(H2s))) & 
                                        indepDatasetsMeta[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+indepWaldLRTCompareMeta[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasetsMeta)[4+indepWaldLRTCompareMeta[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.power) <- c(colnames(indepDatasetsMeta)[1:4], "power", "Type", "hSquare")
  result.detail.power <- as.data.frame(result.detail.power)
  result.detail.power$fdr <- as.numeric(as.vector(result.detail.power$fdr))
  result.detail.power$power <- as.numeric(as.vector(result.detail.power$power))
  
  # a matrix to store fdp estimates
  result.detail.fdp <- matrix(NA, nrow=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)*length(H2s), ncol=7)
  for (h2Indx in 1:length(H2s)) {
    result.detail.fdp[((h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)+1):(h2Indx*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta)), 7] <- rep(H2s[h2Indx], times=length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta))
    for (replicateIndx in 1:length(indepWaldLRTCompareMeta)) {
      for (fdrLvl in 1:length(seq(0,0.2,by=0.01))) {
        result.detail.fdp[(h2Indx-1)*length(seq(0,0.2,by=0.01))*length(indepWaldLRTCompareMeta) + (replicateIndx-1)*length(seq(0,0.2,by=0.01)) + fdrLvl,  1:6] <- 
          unlist(c(indepDatasetsMeta[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasetsMeta)/length(H2s)), 
                                       rep(TRUE, times=nrow(indepDatasetsMeta)/length(H2s)), 
                                       rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasetsMeta)/length(H2s))) & 
                                       indepDatasetsMeta[ , "fdr"]==seq(0,0.2,by=0.01)[fdrLvl], 1:4][1, ], 
                   mean(indepDatasetsMeta[c(rep(FALSE, times=(h2Indx-1)*nrow(indepDatasetsMeta)/length(H2s)), 
                                            rep(TRUE, times=nrow(indepDatasetsMeta)/length(H2s)), 
                                            rep(FALSE, times=(length(H2s)-h2Indx)*nrow(indepDatasetsMeta)/length(H2s))) & 
                                            indepDatasetsMeta[ , "fdr"] == seq(0,0.2,by=0.01)[fdrLvl], 4+(ncol(indepDatasetsMeta)-4)/2+indepWaldLRTCompareMeta[replicateIndx]], na.rm=TRUE), 
                   colnames(indepDatasetsMeta)[4+(ncol(indepDatasetsMeta)-4)/2+indepWaldLRTCompareMeta[replicateIndx]]))
      }
    }
  }
  colnames(result.detail.fdp) <- c(colnames(indepDatasetsMeta)[1:4], "fdp", "Type", "hSquare")
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
  
  b <- ggplot(combinedData[combinedData$Type %in% c("GhostKnockoff, score test","GhostKnockoff, meta score test (Fisher's method)"), ], aes(x=as.numeric(fdr), y=as.numeric(value), group=Type)) + 
    scale_x_continuous("False discovery rate", breaks=seq(from=0, to=0.2, by=0.02)) +
    geom_line(aes(color=Type), size=2) +
    facet_grid_sc(rows=vars(valueType), scales=list(y=scales_y), switch="y", labeller=as_labeller(c(FDP="Observed FDR", Power="Power"))) + 
    geom_abline(data=data.frame(yint=0,yslope=1,valueType="FDP"), aes(intercept=yint, slope=yslope), col='grey', size=2) +
    #ggtitle("Quantitative phenotypes") + 
    ggtitle("Dichotomous phenotypes") +
    scale_color_manual(values=c("#CC79A7", "#fe9929")) +
    labs(color="Methods") +
    ylab(NULL) +
    theme_bw() +    
    theme(legend.position = c(0.37, 0.87),
          plot.title = element_text(face = "bold", size = 12),
          axis.ticks = element_line(colour = "grey70", size = 0.2),
          panel.grid.major = element_line(colour = "grey70", size = 0.2),
          panel.grid.minor = element_blank())

  return(list(indepCombinePlots=indepCombinePlots,
              indepPlots=indepPlots,
              barCombinePlots=barCombinePlots,
              combinePlots=combinePlots,
              metaPlot=b))
}

 
continuous1 <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                          indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                          indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.1)
SchemeA1 <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.1, Kest=0.383)
SchemeC1 <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.1, Kest=1)
SchemeB1 <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.1, Kest=0.808)

continuous <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                          indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                          indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.2)
SchemeA <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.2, Kest=0.383)
SchemeC <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.2, Kest=1)
SchemeB <- compareH1H0(Thetas=c(0,4,7),EpsSDs=c(8,4,1),Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail),methodsForCompare=c(13,15,19,20),
                       indepDatasets=as.data.frame(indep0.8$singleResult.detail),indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail),
                       indepmethodsForCompare=c(13,16),indepWaldLRTCompare=c(13,16,18,19),indepWaldLRTCompareMeta=c(1,2),TypeOfTest="single variant",H2s=2.5,PCausal=10,FDPSpecify=0.2, Kest=0.808)

sixPanel10 <- ggarrange(continuous1$barCombinePlots, SchemeC1$barCombinePlots,
                        nrow=1, ncol=2, common.legend=TRUE, legend="bottom")
sixPanel10 <- annotate_figure(sixPanel10, top= text_grob(paste0("Robustness in different levels of sample relatedness"), face="bold"))

sixPanel10 <- ggarrange(continuous1$barCombinePlots, SchemeA1$barCombinePlots, 
                        SchemeB1$barCombinePlots, SchemeC1$barCombinePlots,
                        nrow=2, ncol=2, common.legend=TRUE, legend="bottom")
sixPanel10 <- ggarrange(continuous$barCombinePlots, SchemeA$barCombinePlots, 
                        SchemeB$barCombinePlots, SchemeC$barCombinePlots,
                        nrow=2, ncol=2, common.legend=TRUE, legend="bottom")

sixPanel02 <- ggarrange(continuous$combinePlots[[1]], continuous$combinePlots[[2]], continuous$combinePlots[[3]], nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
sixPanel02 <- annotate_figure(sixPanel02, top= text_grob(paste0("Related quantitative phenotypes"), face="bold"))
sixPanel03 <- ggarrange(ggarrange(continuous$indepPlots,blank,nrow=1,ncol=2, legend="right", widths=c(1.28,1)), sixPanel02, nrow=2, ncol=1, heights=c(1,1.15))

sixPanel04 <- ggarrange(SchemeA$combinePlots[[1]], SchemeA$combinePlots[[2]], SchemeA$combinePlots[[3]], nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
sixPanel04 <- annotate_figure(sixPanel04, top= text_grob(paste0("Related dichotomous phenotypes; sampling scheme A"), face="bold"))
sixPanel05 <- ggarrange(ggarrange(SchemeA$indepPlots,blank,nrow=1,ncol=2, legend="right", widths=c(1.28,1)), sixPanel04, nrow=2, ncol=1, heights=c(1,1.15))


sixPanel06 <- ggarrange(SchemeB$combinePlots[[1]], SchemeB$combinePlots[[2]], SchemeB$combinePlots[[3]], nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
sixPanel06 <- annotate_figure(sixPanel06, top= text_grob(paste0("Related dichotomous phenotypes; sampling scheme B"), face="bold"))
sixPanel07 <- ggarrange(SchemeC$combinePlots[[1]], SchemeC$combinePlots[[2]], SchemeC$combinePlots[[3]], nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
sixPanel07 <- annotate_figure(sixPanel07, top= text_grob(paste0("Related dichotomous phenotypes; sampling scheme C"), face="bold"))

sixPanel091 <- ggarrange(sixPanel06, sixPanel07, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")

#meta
sixPanel06 <- ggarrange(continuous1$indepCombinePlots, SchemeC1$indepCombinePlots, nrow=1, ncol=2, common.legend=TRUE, legend="bottom")
sixPanel06 <- annotate_figure(sixPanel06, top= text_grob(paste0("Robustness in different types of association tests"), face="bold"))
sixPanel07 <- ggarrange(continuous1$metaPlot, SchemeC1$metaPlot, nrow=1, ncol=2, common.legend=TRUE, legend="bottom")
sixPanel07 <- annotate_figure(sixPanel07, top= text_grob(paste0("Robustness in meta-analysis of multiple studies"), face="bold"))
sixPanel091 <- ggarrange(sixPanel07, sixPanel06, nrow=2, ncol=1, common.legend=TRUE, legend="bottom", labels="AUTO")


Thetas=c(0,4,7)
EpsSDs=c(8,4,1)
Datasets=rbind.data.frame(result02$singleResult.detail,result42$singleResult.detail,result82$singleResult.detail)
methodsForCompare=c(13,15,19,20)
indepDatasets=as.data.frame(indep0.8$singleResult.detail)
indepDatasetsMeta=as.data.frame(indep0.8meta$singleResult.detail)
indepmethodsForCompare=c(13,16)
indepWaldLRTCompare=c(13,16,18,19)#for different types of marginal association tests
indepWaldLRTCompareMeta=c(1,2)
TypeOfTest="single variant"
H2s=2.5
PCausal=10
FDPSpecify=0.1
Kest <- 1
