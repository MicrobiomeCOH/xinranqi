# functions of simulating related genotypes from 10 family members pedigree
getCommonSnipsFromLegendFiles <- function(legendFiles){
  snips <- list(0)
  for (i in 1:length(legendFiles)) {
    snips[[i]] <- sort(read.table(legendFiles[[i]], header=TRUE)[,2])
  }
  snipsTogether <- snips[[1]]
  for (i in 2:length(legendFiles)) {
    snipsTogether <- c(snipsTogether, snips[[i]])
  }
  commonSnipPositionsInChroms <- snipsTogether
  for (i in 1:(length(legendFiles)-1)) {
    commonSnipPositionsInChroms <- commonSnipPositionsInChroms[duplicated(commonSnipPositionsInChroms)]
  }
  commonSnipIndiceInPhasedFiles <- matrix(0,length(legendFiles),length(commonSnipPositionsInChroms))
  for(i in 1:length(legendFiles)){
    commonSnipIndiceInPhasedFiles[i,] <- order(order(c(snips[[i]], commonSnipPositionsInChroms)))[(length(snips[[i]])+1):(length(snips[[i]])+length(commonSnipPositionsInChroms))]-1:length(commonSnipPositionsInChroms)
  }
  list(commonSnipPositionsInChroms=commonSnipPositionsInChroms, commonSnipIndiceInPhasedFiles=commonSnipIndiceInPhasedFiles)
}


generateIndependent <- function(){
  numIndependentAll <- numIndependentDiseased+numIndependentControl       ##dependent individuals
  
  a <- matrix(0,nrow=numIndependentAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  disease <- rep(0,numIndependentAll)
  fromWhichPop <- matrix(0,nrow=numIndependentAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  fromWhichChrom <- matrix(0,nrow=numIndependentAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  
  indexIndependentDiseased <- 0
  indexIndependentControl <- 0
  while ((indexIndependentDiseased < numIndependentDiseased) | (indexIndependentControl < numIndependentControl)) {
    temp.list <- generateIndependentSingle(selectedSnipPositionRange,adjustedRecombinationRate,heritageWeights,selectedSnipsFromPhasedFiles,selectedSnipPositionsInChroms)
    temp1.list <- generateIndependentSingle(selectedSnipPositionRange,adjustedRecombinationRate,heritageWeights,selectedSnipsFromPhasedFiles,selectedSnipPositionsInChroms)
    if (isHNull) {
      heritages <- (apply(tempPopMatrix == temp.list$fromWhichPop, 1, "sum")+apply(tempPopMatrix == temp1.list$fromWhichPop, 1, "sum"))/(2*numOfSelectedSnips)
      thisDiseased <- rbinom(1,1, sum(heritages * H0DiseaseRates))
    } else {
      thisDiseased <- rbinom(1,1,thisDiseaseRate[1+temp.list$result[diseaseSnip]+temp1.list$result[diseaseSnip]])
    }
    if (thisDiseased & (indexIndependentDiseased < numIndependentDiseased)){
      indexIndependentDiseased <- indexIndependentDiseased + 1
      a[2*numIndependentControl+(indexIndependentDiseased*2-1):(indexIndependentDiseased*2),] <- rbind(temp.list$result,temp1.list$result)
      fromWhichPop[2*numIndependentControl+(indexIndependentDiseased*2-1):(indexIndependentDiseased*2),] <- rbind(temp.list$fromWhichPop,temp1.list$fromWhichPop)
      fromWhichChrom[2*numIndependentControl+(indexIndependentDiseased*2-1):(indexIndependentDiseased*2),] <- rbind(temp.list$fromWhichChrom,temp1.list$fromWhichChrom)
      disease[numIndependentControl+indexIndependentDiseased] <- 1
    } else if (!thisDiseased & (indexIndependentControl < numIndependentControl)){
      indexIndependentControl <- indexIndependentControl + 1
      a[(indexIndependentControl*2-1):(indexIndependentControl*2),] <- rbind(temp.list$result,temp1.list$result)
      fromWhichPop[(indexIndependentControl*2-1):(indexIndependentControl*2),] <- rbind(temp.list$fromWhichPop,temp1.list$fromWhichPop)
      fromWhichChrom[(indexIndependentControl*2-1):(indexIndependentControl*2),] <- rbind(temp.list$fromWhichChrom,temp1.list$fromWhichChrom)
      disease[indexIndependentControl] <- 0
    }
  }
  list(a=a,fromWhichPop=fromWhichPop,fromWhichChrom=fromWhichChrom,disease=disease)
}


##This function generates one pair of child chromosomes from his parents (each has one pair chromosomes).
##inputs: 
##positionRange
##rate: recombination rate
##parents: matrix of parent chromosomes, 4 rows, first two are father chromosome, last two mothers.
##chosenSnips: M columns, one for each snip that we are interested in chosenSnips
##output: matrix of child chromosomes, 2 rows, one for each chromosomes, M columns, one for each snip 
generateChild <- function(positionRange,rate,parents.list,chosenSnips){
  parentsSnips <- rbind(parents.list[[1]]$result, parents.list[[2]]$result, parents.list[[3]]$result, parents.list[[4]]$result)
  M <- length(parents.list[[1]]$result)
  result <- matrix(rep(0,2*M),nrow=2,ncol=M)
  fromWhichParent <- matrix(rep(0,2*M),nrow=2,ncol=M) 
  fromWhichPop <- matrix(rep(0,2*M),nrow=2,ncol=M) 
  fromWhichChrom <- matrix(rep(0,2*M),nrow=2,ncol=M) 
  ##generate the chromosome from father
  recombinationNumber <- rpois(1,rate)   #number of recombinations
  ##find end points of those sections (number of end points = number of sections + 1)
  if (recombinationNumber != 0){
    recoms <- sort(floor(runif(recombinationNumber)*(positionRange[2]-positionRange[1]))+min(chosenSnips))
    sectionEnds <- unique(c(1,(order(order(c(chosenSnips,recoms)))[(M+1):(M+recombinationNumber)]-0:(recombinationNumber-1)),M+1))
  } else {
    sectionEnds <- c(1,M+1)
  }
  sections <- rbinom(length(sectionEnds)-1,1,0.5)
  for (j in 1:length(sections)){
    result[1,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parentsSnips[sections[j]+1,sectionEnds[j]:(sectionEnds[j+1]-1)]
    fromWhichParent[1,sectionEnds[j]:(sectionEnds[j+1]-1)] <- sections[j]+1
    fromWhichPop[1,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parents.list[[sections[j]+1]]$fromWhichPop[sectionEnds[j]:(sectionEnds[j+1]-1)]
    fromWhichChrom[1,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parents.list[[sections[j]+1]]$fromWhichChrom[sectionEnds[j]:(sectionEnds[j+1]-1)]
  }
  #chromosome from mother
  recombinationNumber <- rpois(1,rate)
  if (recombinationNumber != 0){
    recoms <- sort(floor(runif(recombinationNumber)*(positionRange[2]-positionRange[1]))+min(chosenSnips))
    sectionEnds <- unique(c(1,(order(order(c(chosenSnips,recoms)))[(M+1):(M+recombinationNumber)]-0:(recombinationNumber-1)),M+1))
  } else {
    sectionEnds <- c(1,M+1)
  }
  sections <- rbinom(length(sectionEnds)-1,1,0.5)
  for (j in 1:length(sections)){
    result[2,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parentsSnips[sections[j]+3,sectionEnds[j]:(sectionEnds[j+1]-1)]
    fromWhichParent[2,sectionEnds[j]:(sectionEnds[j+1]-1)] <- sections[j]+3
    fromWhichPop[2,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parents.list[[sections[j]+3]]$fromWhichPop[sectionEnds[j]:(sectionEnds[j+1]-1)]
    fromWhichChrom[2,sectionEnds[j]:(sectionEnds[j+1]-1)] <- parents.list[[sections[j]+3]]$fromWhichChrom[sectionEnds[j]:(sectionEnds[j+1]-1)]
  }
  list(result=result,fromWhichParent=fromWhichParent,fromWhichPop=fromWhichPop,fromWhichChrom=fromWhichChrom)
}


generateIndependentSingle <- function(selectedSnipPositionRange,adjustedRecombinationRate,heritageWeights,selectedSnipsFromPhasedFiles,selectedSnipPositionsInChroms){
  M <- length(selectedSnipPositionsInChroms)
  result <- rep(0,M)
  fromWhichPop <- rep(0,M)
  fromWhichChrom <- rep(0,M)
  recombinationNumber <- rpois(1,adjustedRecombinationRate)
  if (recombinationNumber != 0){
    recoms <- sort(floor(runif(recombinationNumber)*(selectedSnipPositionRange[2]-selectedSnipPositionRange[1]))+min(selectedSnipPositionsInChroms))
    sectionEnds <- unique(c(1,(order(order(c(selectedSnipPositionsInChroms,recoms)))[(M+1):(M+recombinationNumber)]-0:(recombinationNumber-1)),M+1))
  } else {
    sectionEnds <- c(1,M+1)
  }
  sections <- sample(length(selectedSnipsFromPhasedFiles),length(sectionEnds)-1,replace=TRUE,heritageWeights)
  for (j in 1:length(sections)){
    which <- sample(length(selectedSnipsFromPhasedFiles[[sections[j]]][,1]),1)
    result[sectionEnds[j]:(sectionEnds[j+1]-1)] <- selectedSnipsFromPhasedFiles[[sections[j]]][which,sectionEnds[j]:(sectionEnds[j+1]-1)]
    fromWhichPop[sectionEnds[j]:(sectionEnds[j+1]-1)] <- sections[j]
    fromWhichChrom[sectionEnds[j]:(sectionEnds[j+1]-1)] <- which 
  }
  list(result=result,fromWhichPop=fromWhichPop,fromWhichChrom=fromWhichChrom)
}


generateFamily <- function(){
  numMembersAll <- sum(numMembersPerFamily)
  a <- matrix(0,nrow=numMembersAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  disease <- rep(0,numMembersAll)
  fromWhichPop <- matrix(0,nrow=numMembersAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  fromWhichChrom <- matrix(0,nrow=numMembersAll*2,ncol=numOfSelectedSnips) ##An empty matrix to take the grand result
  
  #result of family members
  indexOfMembers <- 0
  for(i in 1:numFamilies){
    nonChildIsDiseased <- TRUE
    while (nonChildIsDiseased) {
      diseaseInFamily <- rep(0,numMembersPerFamily[i])
      parents.list <- list()
      #generate 4 independent chromosomes, two for each parents
      for (j in 1:4){
        parents.list <- c(parents.list, list(generateIndependentSingle(selectedSnipPositionRange,adjustedRecombinationRate,heritageWeights,selectedSnipsFromPhasedFiles,selectedSnipPositionsInChroms)))
      }
      #generate children of the two parents
      tempC.list <- list()
      for(j in 1:2){
        tempC.list <- c(tempC.list, list(generateChild(selectedSnipPositionRange,0,parents.list,selectedSnipPositionsInChroms)))
      }
      # second generation 4 independent chromosome, corresponding to above children
      secondParents.list <- list()
      #generate 4 independent chromosomes, two for each parents
      for (j in 1:4){
        secondParents.list <- c(secondParents.list, list(generateIndependentSingle(selectedSnipPositionRange,adjustedRecombinationRate,heritageWeights,selectedSnipsFromPhasedFiles,selectedSnipPositionsInChroms)))
      }
      # third generation grand-children
      tempC2.list <- list()
      for (j in 1:2) {
        tempP.list <- list()
        tempP.list <- c(tempP.list, secondParents.list[(j-1)*2+1], secondParents.list[j*2])
        tempP.list[[3]] <- list(result=tempC.list[[j]]$result[1,], fromWhichPop=tempC.list[[j]]$fromWhichPop[1,], fromWhichChrom=tempC.list[[j]]$fromWhichChrom[1,])
        tempP.list[[4]] <- list(result=tempC.list[[j]]$result[2,], fromWhichPop=tempC.list[[j]]$fromWhichPop[2,], fromWhichChrom=tempC.list[[j]]$fromWhichChrom[2,])
        for (jThird in 1:2) {
          tempC2.list <- c(tempC2.list, list(generateChild(selectedSnipPositionRange,0,tempP.list,selectedSnipPositionsInChroms)))
        }
      }
      
      # fulfillment of genotypes, disease status, population & chromosomes
      for (j in 1:4){
        # grandparents
        a[indexOfMembers*2+j,] <- parents.list[[j]]$result
        fromWhichPop[indexOfMembers*2+j,] <- parents.list[[j]]$fromWhichPop
        fromWhichChrom[indexOfMembers*2+j,] <- parents.list[[j]]$fromWhichChrom
        # independent second generation parents
        a[indexOfMembers*2+8+j,] <- secondParents.list[[j]]$result
        fromWhichPop[indexOfMembers*2+8+j,] <- secondParents.list[[j]]$fromWhichPop
        fromWhichChrom[indexOfMembers*2+8+j,] <- secondParents.list[[j]]$fromWhichChrom
        # grandchildren
        a[indexOfMembers*2+12+(j*2-1):(j*2),] <- tempC2.list[[j]]$result
        fromWhichPop[indexOfMembers*2+12+(j*2-1):(j*2),] <- tempC2.list[[j]]$fromWhichPop
        fromWhichChrom[indexOfMembers*2+12+(j*2-1):(j*2),] <- tempC2.list[[j]]$fromWhichChrom
      }
      for (j in 1:2){
        # second generation parents
        a[indexOfMembers*2+4+(j*2-1):(j*2),] <- tempC.list[[j]]$result
        fromWhichPop[indexOfMembers*2+4+(j*2-1):(j*2),] <- tempC.list[[j]]$fromWhichPop
        fromWhichChrom[indexOfMembers*2+4+(j*2-1):(j*2),] <- tempC.list[[j]]$fromWhichChrom
      }
      # hyper-parameters
      for(j in 1:numMembersPerFamily[i]){
        if (isHNull) {
          heritages <- (apply(tempPopMatrix == fromWhichPop[indexOfMembers*2+2*j-1,], 1, "sum")+apply(tempPopMatrix == fromWhichPop[indexOfMembers*2+2*j,], 1, "sum"))/(2*numOfSelectedSnips)
          thisDiseased <- rbinom(1,1, sum(heritages * H0DiseaseRates))
        } else {
          thisDiseased <- rbinom(1,1,thisDiseaseRate[1+a[indexOfMembers*2+2*j-1,diseaseSnip]+a[indexOfMembers*2+2*j,diseaseSnip]])
        }
        diseaseInFamily[j] <- thisDiseased
      }
      nonChildIsDiseased <- sum(diseaseInFamily[-c(1:2, 5:6)]) == 0
    }
    disease[indexOfMembers+1:numMembersPerFamily[i]] <- diseaseInFamily
    indexOfMembers <- indexOfMembers+numMembersPerFamily[i]
  }
  list(a=a,fromWhichPop=fromWhichPop,fromWhichChrom=fromWhichChrom,disease=disease)
}


# convert R matrix to Latex matrix
bmatrix <- function(x, digits=NULL, ...) {
  library(xtable)
  default_args = list(include.colnames=FALSE, only.contents=TRUE,
                      include.rownames=FALSE, hline.after=NULL, comment=FALSE,
                      print.results=FALSE)
  passed_args = list(...)
  calling_args = c(list(x=xtable(x, digits=digits)),
                   c(passed_args,
                     default_args[setdiff(names(default_args), names(passed_args))]))
  cat("\\begin{bmatrix}\n",
      do.call(print.xtable, calling_args),
      "\\end{bmatrix}\n")
}


# a function to calculate GLM p values
Get.p.modifiedGLM <- function(G, obj_nullmodel_glm){
  G <- as.matrix(G)
  mu <- obj_nullmodel_glm$fitted.values
  Y.res <- obj_nullmodel_glm$y - mu
  # mu<-result.prelim$nullglm$fitted.values; Y.res<-result.prelim$Y-mu
  outcome <-ifelse(obj_nullmodel_glm$family[1]=="binomial", "D", "C")
  # outcome<-result.prelim$out_type
  if(outcome=='D'){
    if (length(obj_nullmodel_glm$model)==1) {
      p <- ScoreTest_SPA(genos=t(G), pheno=obj_nullmodel_glm$y, cov=NULL, method="fastSPA", minmac=-Inf)$p.value
    } else {
      p <- ScoreTest_SPA(t(G), obj_nullmodel_glm$y, obj_nullmodel_glm$model["X1"], method=c("fastSPA"), minmac=-Inf)$p.value
    }
    #p <- ScoreTest_SPA(t(X),result.prelim$Y,result.prelim$X,method=c("fastSPA"),minmac=-Inf)$p.value
  } else {
    v <- rep(as.numeric(var(Y.res)), nrow(G))
    # required parameters calculated from KS.prelim() function
    if (length(obj_nullmodel_glm$model)==1) {
      X0 <- NULL
    } else {
      X0 <- svd(as.matrix(obj_nullmodel_glm$model["X1"]))$u
    }
    X0 <- cbind(rep(1,length(obj_nullmodel_glm$y)),X0)
    # prepare inverse matrix for covariates
    inv.X0 <- solve(t(X0)%*%(v*X0))
    p <- pchisq((t(G)%*%Y.res)^2/(apply(G*(v*G),2,sum)-apply(t(G)%*%(v*X0)%*%inv.X0*t(t(X0)%*%as.matrix(v*G)),1,sum)),df=1,lower.tail=F)
    # p<-pchisq((t(X)%*%Y.res)^2/(apply(X*(v*X),2,sum)-apply(t(X)%*%(v*result.prelim$X0)%*%result.prelim$inv.X0*t(t(result.prelim$X0)%*%as.matrix(v*X)),1,sum)),df=1,lower.tail=F)
  }
  p[is.na(p)] <- 1
  return(as.matrix(p))
}


# modify STAAR function to include negative MAF for knockoff counterparts
STAAR.modified <- function(genotype,obj_nullmodel,annotation_phred=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,yResponse,out_type){
  
  if(class(genotype) != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }
  
  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }
  
  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & dim(genotype)[2] != dim(annotation_phred)[1]){
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }
  
  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector(MAF<rare_maf_cutoff)
  #RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]
  
  rm(genotype)
  gc()
  annotation_phred <- annotation_phred[RV_label,,drop=FALSE]
  
  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    MAF <- MAF[RV_label]
    rm(Geno_rare)
    gc()
    
    annotation_rank <- 1 - 10^(-annotation_phred/10)
    
    ## beta(1,25)
    w_1 <- dbeta(MAF,1,25)
    ## beta(1,1)
    w_2 <- dbeta(MAF,1,1)
    if(dim(annotation_phred)[2] == 0){
      ## Burden, SKAT, ACAT-V
      w_B <- w_S <- as.matrix(cbind(w_1,w_2))
      #w_A <- as.matrix(cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_2^2/dbeta(MAF,0.5,0.5)^2))
      w_A <- as.matrix(cbind(sapply(1:length(MAF), function(x) {ifelse(MAF[x]<=0, 0, w_1[x]^2 / dbeta(MAF[x],0.5,0.5)^2)}),
                             sapply(1:length(MAF), function(x) {ifelse(MAF[x]<=0, 0, w_2[x]^2 / dbeta(MAF[x],0.5,0.5)^2)})))
    }else{
      ## Burden
      w_B_1 <- annotation_rank*w_1
      w_B_1 <- cbind(w_1,w_B_1)
      w_B_2 <- annotation_rank*w_2
      w_B_2 <- cbind(w_2,w_B_2)
      w_B <- cbind(w_B_1,w_B_2)
      w_B <- as.matrix(w_B)
      
      ## SKAT
      w_S_1 <- sqrt(annotation_rank)*w_1
      w_S_1 <- cbind(w_1,w_S_1)
      w_S_2 <- sqrt(annotation_rank)*w_2
      w_S_2 <- cbind(w_2,w_S_2)
      w_S <- cbind(w_S_1,w_S_2)
      w_S <- as.matrix(w_S)
      
      ## ACAT-V
      w_A_1 <- annotation_rank*w_1^2/dbeta(MAF,0.5,0.5)^2
      w_A_1 <- cbind(w_1^2/dbeta(MAF,0.5,0.5)^2,w_A_1)
      w_A_2 <- annotation_rank*w_2^2/dbeta(MAF,0.5,0.5)^2
      w_A_2 <- cbind(w_2^2/dbeta(MAF,0.5,0.5)^2,w_A_2)
      w_A <- cbind(w_A_1,w_A_2)
      w_A <- as.matrix(w_A)
    }
    
    if(obj_nullmodel$relatedness){
      if(!obj_nullmodel$sparse_kins){
        P <- obj_nullmodel$P
        P_scalar <- sqrt(dim(P)[1])
        P <- P*P_scalar
        
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)
        Mac <- as.integer(round(MAF*2*dim(G)[1]))
        Mac[Mac<0] <- 0
        
        pvalues <- STAAR::STAAR_O_SMMAT(G,P,residuals.phenotype,
                                        weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                        mac=Mac, mac_thres=-1)
      }else{
        Sigma_i <- obj_nullmodel$Sigma_i
        Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
        cov <- obj_nullmodel$cov
        
        residuals.phenotype <- obj_nullmodel$scaled.residuals
        Mac <- as.integer(round(MAF*2*dim(G)[1]))
        Mac[Mac<0] <- 0
        
        pvalues <- STAAR::STAAR_O_SMMAT_sparse(G,Sigma_i,Sigma_iX,cov,residuals.phenotype,
                                               weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                               mac=Mac, mac_thres=-1)
      }
    }else{
      X <- model.matrix(obj_nullmodel)
      working <- obj_nullmodel$weights
      sigma <- sqrt(summary(obj_nullmodel)$dispersion)
      if(obj_nullmodel$family[1] == "binomial"){
        fam <- 1
      }else if(obj_nullmodel$family[1] == "gaussian"){
        fam <- 0
      }
      
      residuals.phenotype <- obj_nullmodel$y - obj_nullmodel$fitted.values
      Mac <- as.integer(round(MAF*2*dim(G)[1]))
      Mac[Mac<0] <- 0
      
      pvalues <- STAAR::STAAR_O(G,X,working,sigma,fam,residuals.phenotype,
                                weights_B=w_B,weights_S=w_S,weights_A=w_A,
                                mac=Mac, mac_thres=-1)
    }
    
    # add SKAT(1,25) and updated SKAT(1,1) if the first element is NaN due to MAF==0.5 elementwise
    if (any(is.na(pvalues[1:2]))) {
      SKATobj <- if (length(obj_nullmodel$coefficients)==1) {SKAT_Null_Model(yResponse~1, out_type=out_type)} else {SKAT_Null_Model(yResponse~as.matrix(obj_nullmodel$X[,-1]), out_type=out_type)}
      for (SKATi in 1:2) {
        pvalues[SKATi] <- SKAT(G, SKATobj, weights.beta=c(1, c(25, 1)[SKATi]), r.corr=0)$p.value
      }
    }
    
    # add Burden(1,25) & Burden(1,1) if fit_null_glmmkin reduce to fit_null_glm and p-values are all NAs
    if (any(is.na(pvalues[3:4]))) {
      Burdenobj <- if (length(obj_nullmodel$coefficients)==1) {SKAT_Null_Model(yResponse~1, out_type=out_type)} else {SKAT_Null_Model(yResponse~as.matrix(obj_nullmodel$X[,-1]), out_type=out_type)}
      for (Burdeni in 1:2) {
        pvalues[2+Burdeni] <- SKAT(G, Burdenobj, weights.beta=c(1, c(25, 1)[Burdeni]), r.corr=1)$p.value
      }
    }
    
    # add ACAT-V(1,25) & ACAT-V(1,11) if fit_null_glmmkin reduce to fit_null_glm and p-values are all NAs
    if (any(is.na(pvalues[5:6]))) {
      ACATobj <- if (length(obj_nullmodel$coefficients)==1) {NULL_Model(as.vector(yResponse))} else {NULL_Model(as.vector(yResponse), Z=as.matrix(obj_nullmodel$X[,-1]))}
      for (ACATi in 1:2) {
        pvalues[4+ACATi] <- ACAT_V(G, ACATobj, weights.beta=c(1, c(25, 1)[ACATi]))
      }
    } 
    
    num_variant <- sum(RV_label) #dim(G)[2]
    cMAC <- sum(G)
    num_annotation <- dim(annotation_phred)[2]+1
    results_STAAR_O <- STAAR::CCT(pvalues)
    results_ACAT_O <- STAAR::CCT(pvalues[c(1,num_annotation+1,2*num_annotation+1,3*num_annotation+1,4*num_annotation+1,5*num_annotation+1)])
    pvalues_STAAR_S_1_25 <- STAAR::CCT(pvalues[1:num_annotation])
    pvalues_STAAR_S_1_1 <- STAAR::CCT(pvalues[(num_annotation+1):(2*num_annotation)])
    pvalues_STAAR_B_1_25 <- STAAR::CCT(pvalues[(2*num_annotation+1):(3*num_annotation)])
    pvalues_STAAR_B_1_1 <- STAAR::CCT(pvalues[(3*num_annotation+1):(4*num_annotation)])
    pvalues_STAAR_A_1_25 <- STAAR::CCT(pvalues[(4*num_annotation+1):(5*num_annotation)])
    pvalues_STAAR_A_1_1 <- STAAR::CCT(pvalues[(5*num_annotation+1):(6*num_annotation)])
    
    results_STAAR_S_1_25 <- c(pvalues[1:num_annotation],pvalues_STAAR_S_1_25)
    results_STAAR_S_1_25 <- data.frame(t(results_STAAR_S_1_25))
    
    results_STAAR_S_1_1 <- c(pvalues[(num_annotation+1):(2*num_annotation)],pvalues_STAAR_S_1_1)
    results_STAAR_S_1_1 <- data.frame(t(results_STAAR_S_1_1))
    
    results_STAAR_B_1_25 <- c(pvalues[(2*num_annotation+1):(3*num_annotation)],pvalues_STAAR_B_1_25)
    results_STAAR_B_1_25 <- data.frame(t(results_STAAR_B_1_25))
    
    results_STAAR_B_1_1 <- c(pvalues[(3*num_annotation+1):(4*num_annotation)],pvalues_STAAR_B_1_1)
    results_STAAR_B_1_1 <- data.frame(t(results_STAAR_B_1_1))
    
    results_STAAR_A_1_25 <- c(pvalues[(4*num_annotation+1):(5*num_annotation)],pvalues_STAAR_A_1_25)
    results_STAAR_A_1_25 <- data.frame(t(results_STAAR_A_1_25))
    
    results_STAAR_A_1_1 <- c(pvalues[(5*num_annotation+1):(6*num_annotation)],pvalues_STAAR_A_1_1)
    results_STAAR_A_1_1 <- data.frame(t(results_STAAR_A_1_1))
    
    if(dim(annotation_phred)[2] == 0){
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)","STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)","STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)","STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)","STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)","STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)","STAAR-A(1,1)")
    }else{
      colnames(results_STAAR_S_1_25) <- c("SKAT(1,25)",
                                          paste0("SKAT(1,25)-",colnames(annotation_phred)),
                                          "STAAR-S(1,25)")
      colnames(results_STAAR_S_1_1) <- c("SKAT(1,1)",
                                         paste0("SKAT(1,1)-",colnames(annotation_phred)),
                                         "STAAR-S(1,1)")
      colnames(results_STAAR_B_1_25) <- c("Burden(1,25)",
                                          paste0("Burden(1,25)-",colnames(annotation_phred)),
                                          "STAAR-B(1,25)")
      colnames(results_STAAR_B_1_1) <- c("Burden(1,1)",
                                         paste0("Burden(1,1)-",colnames(annotation_phred)),
                                         "STAAR-B(1,1)")
      colnames(results_STAAR_A_1_25) <- c("ACAT-V(1,25)",
                                          paste0("ACAT-V(1,25)-",colnames(annotation_phred)),
                                          "STAAR-A(1,25)")
      colnames(results_STAAR_A_1_1) <- c("ACAT-V(1,1)",
                                         paste0("ACAT-V(1,1)-",colnames(annotation_phred)),
                                         "STAAR-A(1,1)")
    }
    
    return(list(num_variant = num_variant,
                cMAC = cMAC,
                RV_label = RV_label,
                results_STAAR_O = results_STAAR_O,
                results_ACAT_O = results_ACAT_O,
                results_STAAR_S_1_25 = results_STAAR_S_1_25,
                results_STAAR_S_1_1 = results_STAAR_S_1_1,
                results_STAAR_B_1_25 = results_STAAR_B_1_25,
                results_STAAR_B_1_1 = results_STAAR_B_1_1,
                results_STAAR_A_1_25 = results_STAAR_A_1_25,
                results_STAAR_A_1_1 = results_STAAR_A_1_1))
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }
}


# a function to estimate inflation factor \eta
etaEstimate <- function(N, ThetaHat, Var.e.Hat, KinshipMatrix) {
  sqrt(N/ (N + 2*ThetaHat*sum(KinshipMatrix[upper.tri(KinshipMatrix, diag=FALSE)]^2)/(ThetaHat+Var.e.Hat) ))
}


ES.prelim.modified <- function(cor.G,sd.G=NULL,M=1,method='asdp',max.size=500, EtaHat){
  if(length(sd.G)!=0){temp.index<-which(sd.G!=0)}else{temp.index<-1:nrow(cor.G)}
  n.G<-nrow(cor.G)
  #Permutation test for constant variants in the reference panel
  permute.index<-rep(0,n.G)
  permute.index[-temp.index]<-1
  
  Normal_50Studies<-matrix(rnorm(n.G*M*50),n.G*M,50)
  P.each<-matrix(0,n.G,n.G)
  if(length(temp.index)!=0){
    Sigma<-cor.G[temp.index,temp.index,drop=F]
    SigmaInv<-solve(Sigma)#invcov.shrink(Sigma,verbose=F)
    
    if(method=='sdp'){temp.s<-create.solve_sdp_M(Sigma,M=M)}
    if(method=='asdp'){temp.s<-create.solve_asdp_M(Sigma,M=M,max.size=max.size)}
    s<-temp.s
    diag_s<-diag(s,length(s))
    
    if(sum(s)==0){
      V.left<-matrix(0,length(temp.index)*M,length(temp.index)*M)
    }else{
      Sigma_k<-2*diag_s - s*t(s*SigmaInv)
      V.each<-Matrix(forceSymmetric(Sigma_k-diag_s))
      
      #random part of knockoff
      V<-matrix(1,M,M)%x%V.each
      diag(V)<-diag(V)+rep(s,M)
      V.left<-try(t(chol(V,pivot = TRUE,tol=1e-6)),silent=T)
      if(class(V.left)=="try-error"){
        svd.fit<-svd(V)
        u<-svd.fit$u
        svd.fit$d[is.na(svd.fit$d)]<-0
        cump<-cumsum(svd.fit$d)/sum(svd.fit$d)
        n.svd<-which(cump>=0.999)[1]
        if(is.na(n.svd)){n.svd<-nrow(V)}
        svd.index<-intersect(1:n.svd,which(svd.fit$d!=0))
        V.left<-t(sqrt(svd.fit$d[svd.index])*t(u[,svd.index,drop=F]))
      }
    }
    P.each[temp.index,temp.index]<-diag(1,length(s))-s*SigmaInv
    V.index<-rep(temp.index,M)+rep(0:(M-1),each=length(temp.index))*n.G
    Normal_50Studies[V.index,]<-as.matrix(EtaHat*(V.left%*%matrix(rnorm(ncol(V.left)*50),ncol(V.left),50)))
    # Normal_50Studies[V.index,]<-as.matrix(V.left%*%matrix(rnorm(ncol(V.left)*50),ncol(V.left),50))
    
    #Permutation test for tightly linked variants
    permute.index[temp.index[s==0]]<-1
  }
  permute.V.index<-rep(permute.index,M)
  P.each[permute.index==1,]<-0
  Normal_50Studies[permute.V.index==1,]<-matrix(rnorm(sum(permute.index)*M*50),sum(permute.index)*M,50)
  
  return(list(P.each=as.matrix(P.each), Normal_50Studies=as.matrix(Normal_50Studies), permute.index=permute.index, M=M))
}
