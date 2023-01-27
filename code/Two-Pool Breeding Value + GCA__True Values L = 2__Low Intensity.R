library(AlphaSimR)
Rcpp::sourceCpp("GPCP_polyploids.cpp")

args = commandArgs(trailingOnly=TRUE)
STARTERS = as.character(args[[1]])



maxAvoidPlan = function(nInd, nProgeny=2L){
  crossPlan = matrix(1:nInd, ncol=2, byrow=TRUE)
  tmp = c(seq(1, nInd, by=2),
          seq(2, nInd, by=2))
  crossPlan = cbind(rep(tmp[crossPlan[,1]], 
                        each=nProgeny),
                    rep(tmp[crossPlan[,2]], 
                        each=nProgeny))
  return(crossPlan)
}



fullInbredPopG <- function(pop, ploidy){  
  pfreq <- colMeans(pullQtlGeno(pop))/ploidy
  qfreq <- 1 - pfreq
  a <- SP$traits[[1]]@addEff
  Mf <- sum((pfreq - qfreq) * a)
  Mf <- Mf + SP$traits[[1]]@intercept
  return(Mf)
}


calcF = function(geno, p, ploidy){
  # Center genotypes by reference allele frequency
  geno = sweep(x=geno, MARGIN=2, STATS=ploidy*p, FUN="-")
  
  # Calculate the denominator using the expected cumulative 
  # genotype variance under HWE (binomial distribution)
  denom = sum(ploidy*p*(1-p))
  
  G = tcrossprod(geno)/denom
  F = (mean(diag(G))-1) / (ploidy-1)
  return(F)
}





#STARTERS <- list.files(getwd())
#STARTERS <- STARTERS[grep("*.RData", STARTERS)]
#H2 <- c(0.0001, 0.25, 0.5, 0.75)
nCycles = 100
scenario = "TBVGCA_TRUE_maxAvoidCtrlGS"


#for(i in 1:length(H2)){
  #for(i in 3){
  #print(paste("H2 =", H2[i]))
  
  
  #for(j in 1:10){
    #for(j in 1){
    
    #load(STARTERS[j])
    load(STARTERS)
	print(paste(colnames(scenParam), scenParam[1,]))
    
    
    #for(REP in 1:10){
      #for(REP in 1){
      #print(paste("REP", REP))
      
      
      
      output = data.frame(cycle = 1:nCycles,        
                          rep = rep(scenParam$rep, nCycles),
                          scenario = rep(scenario, nCycles),
                          meanDD = rep(scenParam$meanDD, nCycles),
                          varDD = rep(scenParam$varDD, nCycles),
                          ploidy = rep(scenParam$ploidy, nCycles),
                          nQTL = rep(scenParam$nQTL, nCycles),
                          nSNP = rep(scenParam$nSNP, nCycles),
                          nChr = rep(scenParam$nChr, nCycles),
                          H2 = rep(1, nCycles),
                          
                          meanG.Hybrids = numeric(nCycles),
                          meanG.Elite = numeric(nCycles),
                          meanG.EliteA = numeric(nCycles),
                          meanG.EliteB = numeric(nCycles),
                          
                          meanA.Hybrids = numeric(nCycles),
                          meanA.Elite = numeric(nCycles),
                          meanA.EliteA = numeric(nCycles),
                          meanA.EliteB = numeric(nCycles),
                          
                          meanD.Hybrids = numeric(nCycles),
                          meanD.Elite = numeric(nCycles),
                          meanD.EliteA = numeric(nCycles),
                          meanD.EliteB = numeric(nCycles),
                          
                          meanGCA.A = numeric(nCycles),
                          meanGCA.B = numeric(nCycles),
                          meanTrueGCA.A = numeric(nCycles),
                          meanTrueGCA.B = numeric(nCycles),
                          
                          meanTrueSCA.Hybrids = numeric(nCycles),
                          
                          acc.Elite = numeric(nCycles),
                          accGCA.EliteA = numeric(nCycles),
                          accGCA.EliteB = numeric(nCycles),
                          accGPCP.Elite = numeric(nCycles),
                          acc.EliteA = numeric(nCycles),
                          acc.EliteB = numeric(nCycles),
                          
                          varG.Hybrids = numeric(nCycles),
                          varG.Elite = numeric(nCycles),
                          varG.EliteA = numeric(nCycles),
                          varG.EliteB = numeric(nCycles),
                          
                          varA.Hybrids = numeric(nCycles),
                          varA.Elite = numeric(nCycles),
                          varA.EliteA = numeric(nCycles),
                          varA.EliteB = numeric(nCycles),
                          
                          varD.Hybrids = numeric(nCycles),
                          varD.Elite = numeric(nCycles),
                          varD.EliteA = numeric(nCycles),
                          varD.EliteB = numeric(nCycles),
                          
                          varGCA.EliteA = numeric(nCycles),
                          varGCA.EliteB = numeric(nCycles),
                          
                          varSCA.Hybrids = numeric(nCycles),
                          
                          meanInbr.Elite = numeric(nCycles),
                          meanInbr.EliteA = numeric(nCycles),
                          meanInbr.EliteB = numeric(nCycles),
                          meanInbr.Hybrids = numeric(nCycles),
                          
                          meanInbrDep.Elite = numeric(nCycles),
                          meanInbrDep.Hybrids = numeric(nCycles),
                          meanInbrDep.EliteA = numeric(nCycles),
                          meanInbrDep.EliteB = numeric(nCycles),
                          
                          meanPanHet = numeric(nCycles),
                          meanInbrMidHet = numeric(nCycles),
                          
                          mu.Elite = numeric(nCycles),
                          mu.Hybrids = numeric(nCycles),
                          mu.EliteA = numeric(nCycles),
                          mu.EliteB = numeric(nCycles),
                          
                          mu_HW.Hybrids = numeric(nCycles),
                          mu_HW.EliteA = numeric(nCycles),
                          mu_HW.EliteB = numeric(nCycles),
                          
                          stringsAsFactors=FALSE)
      
      
      
      for(cycle in 1:nCycles){
        #for(cycle in 1:10){
        print(paste("Cycle", cycle))
        
        # make the within-pool crosses
        if(cycle == 1){
          EliteA <- newPop(founderPop[1:20])
          EliteB <- newPop(founderPop[21:40])
		      refPop <- c(EliteA, EliteB)
		      pRef = colMeans(pullQtlGeno(refPop))/scenParam$ploidy
          
          EliteA <- randCross(EliteA, nCrosses = 50, nProgeny = 5)
          EliteB <- randCross(EliteB, nCrosses = 50, nProgeny = 5)
		  
		      EliteA <- selectWithinFam(EliteA, nInd = 4, use = "bv")
		      EliteB <- selectWithinFam(EliteB, nInd = 4, use = "bv")
        }
        
        if(cycle > 1){
          
          EliteA <- selectWithinFam(EliteA, nInd = 2, use = "ebv") # select on GCA
          EliteB <- selectWithinFam(EliteB, nInd = 2, use = "ebv") # select on GCA
          
          EliteA <- makeCross(EliteA, crossPlan = maxAvoidPlan(nInd = 100, nProgeny = 5))
          EliteB <- makeCross(EliteB, crossPlan = maxAvoidPlan(nInd = 100, nProgeny = 5))

          # Select within pool on breeding value
          EliteA <- selectWithinFam(EliteA, nInd = 4, use = "bv")
          EliteB <- selectWithinFam(EliteB, nInd = 4, use = "bv")
        }
        
        

        
        
        #half diallel = 1000*1001/2 = 500500...
        # make testcrosses
        TestersA <- selectInd(EliteA, nInd = 2, use = "rand")
        TestersB <- selectInd(EliteB, nInd = 2, use = "rand")
        AxB <- makeCross2(EliteA, TestersB, crossPlan = expand.grid(1:nInd(EliteA), 1:nInd(TestersB)), nProgeny = 1)
        BxA <- makeCross2(TestersA, EliteB, crossPlan = expand.grid(1:nInd(TestersA), 1:nInd(EliteB)), nProgeny = 1)
        
        
        # phenotype the testcrosses
        AxB <- setPheno(AxB, H2 = 1, fixEff = cycle)
        BxA <- setPheno(BxA, H2 = 1, fixEff = cycle)
        Hybrids <- c(AxB, BxA)
        
        
        # use GPCP_polyploids.cpp (CalcCrossMeans.cpp) for ploidy-agnostic GCA calculations
        # could be improved by changing source to make only inter-pool crosses
        crossMeans <- calcCrossMean(pullQtlGeno(c(EliteA, EliteB)),
                       SP$traits[[1]]@addEff,
                       SP$traits[[1]]@domEff,
                       scenParam$ploidy)
        
        # subset the inter-pool crosses only
        crossMeans <- crossMeans[crossMeans[ ,1] %in% c(1:nInd(EliteA)) & 
                                   crossMeans[,2] %in% c((nInd(EliteA) + 1):(nInd(EliteA) + nInd(EliteB))), ]
        
        # put inter-pool cross predicted values in matrix
        # rows are EliteA Parent
        # columns are EliteB Parent
        crossMatrix <- matrix(data = crossMeans[ ,3], 
                              nrow = nInd(EliteA), 
                              ncol = nInd(EliteB), 
                              byrow = TRUE)
        
        # each of these means are ~0
        EliteA.TrueGCA <- rowMeans(crossMatrix - mean(crossMatrix))
        EliteB.TrueGCA <- colMeans(crossMatrix - mean(crossMatrix))
        Hybrids.TrueSCA <- crossMatrix - 
                            mean(crossMatrix) -
                            matrix(data = EliteA.TrueGCA,
                                   nrow = nInd(EliteA),
                                   ncol = nInd(EliteB),
                                   byrow = TRUE) -
                            matrix(data = EliteB.TrueGCA,
                                   nrow = nInd(EliteA),
                                   ncol = nInd(EliteB),
                                   byrow = FALSE)
        
        EliteA@ebv <- as.matrix(EliteA.TrueGCA)
        EliteB@ebv <- as.matrix(EliteB.TrueGCA)
          
 
        
        
        # make the genParam
        GP.Hybrids <- genParam(Hybrids)
        GP.EliteA <- genParam(EliteA)
        GP.EliteB <- genParam(EliteB)
        
        
        # record the mean genetic values
        output[cycle, "meanG.Hybrids"] <- meanG(Hybrids)
        output[cycle, "meanG.EliteA"] <- meanG(EliteA)
        output[cycle, "meanG.EliteB"] <- meanG(EliteB)
        
        # record the mean additive genetic values
        output[cycle, "meanA.Hybrids"] <- mean(GP.Hybrids$gv_a)
        output[cycle, "meanA.EliteA"] <- mean(GP.EliteA$gv_a)
        output[cycle, "meanA.EliteB"] <- mean(GP.EliteB$gv_a)
        
        # record the mean dominance genetic value
        output[cycle, "meanD.Hybrids"] <- mean(GP.Hybrids$gv_d)
        output[cycle, "meanD.EliteA"] <- mean(GP.EliteA$gv_d)
        output[cycle, "meanD.EliteB"] <- mean(GP.EliteB$gv_d)
        
        # record the mean true and estimated GCA
        output[cycle, "meanGCA.A"] <- mean(EliteA@ebv)
        output[cycle, "meanGCA.B"] <- mean(EliteB@ebv)
        output[cycle, "meanTrueGCA.A"] <- mean(EliteA.TrueGCA)
        output[cycle, "meanTrueGCA.B"] <- mean(EliteB.TrueGCA)
        
        # record mean true SCA (of the hybrid combinations)
        output[cycle, "meanTrueSCA.Hybrids"] <- mean(Hybrids.TrueSCA)
        
        
        # record the selection accuracy
        output[cycle, "accGCA.EliteA"] <- cor(EliteA@ebv, EliteA.TrueGCA) #EliteA
        output[cycle, "accGCA.EliteB"] <- cor(EliteB@ebv, EliteB.TrueGCA) #EliteB
        
        # record genetic variances: inbred and hybrid
        output[cycle, "varG.Hybrids"] <- varG(Hybrids)
        output[cycle, "varG.EliteA"] <- varG(EliteA)
        output[cycle, "varG.EliteB"] <- varG(EliteB)
        
        output[cycle, "varA.Hybrids"] <- varA(Hybrids)
        output[cycle, "varA.EliteA"] <- varA(EliteA)
        output[cycle, "varA.EliteB"] <- varA(EliteB)
        
        output[cycle, "varD.Hybrids"] <- varD(Hybrids)
        output[cycle, "varD.EliteA"] <- varD(EliteA)
        output[cycle, "varD.EliteB"] <- varD(EliteB)
        
        
        #record true GCA variance
        output[cycle, "varGCA.EliteA"] <- var(EliteA.TrueGCA)
        output[cycle, "varGCA.EliteB"] <- var(EliteB.TrueGCA)
        
        #record true SCA variance
        output[cycle, "varSCA.Hybrids"] <- var(as.vector(Hybrids.TrueSCA))
        
        # record inbreeding: inbred and hybrid
						output[cycle, "meanInbr.Hybrids"] <- calcF(pullQtlGeno(Hybrids), p = pRef, ploidy = scenParam$ploidy)
		output[cycle, "meanInbr.EliteA"] <- calcF(pullQtlGeno(EliteA), p = pRef, ploidy = scenParam$ploidy)
		output[cycle, "meanInbr.EliteB"] <- calcF(pullQtlGeno(EliteB), p = pRef, ploidy = scenParam$ploidy)
        
        
        # record the inbreeding depression
        # option: make a fake population in HWE and look at change in mean
        output[cycle, "meanInbrDep.Hybrids"] <- fullInbredPopG(Hybrids, scenParam$ploidy)
        output[cycle, "meanInbrDep.EliteA"] <- fullInbredPopG(EliteA, scenParam$ploidy)
        output[cycle, "meanInbrDep.EliteB"] <- fullInbredPopG(EliteB, scenParam$ploidy)
        
        
        # record population mean
        output[cycle, "mu.Hybrids"] <- GP.Hybrids$mu
        output[cycle, "mu.EliteA"] <- GP.EliteA$mu
        output[cycle, "mu.EliteB"] <- GP.EliteB$mu
        
        
        # record HWE population mean
        output[cycle, "mu_HW.Hybrids"] <- GP.Hybrids$mu_HW
        output[cycle, "mu_HW.EliteA"] <- GP.EliteA$mu_HW
        output[cycle, "mu_HW.EliteB"] <- GP.EliteB$mu_HW
        
        
        # record panmictic heterosis
        #mean of hybrids - mean of two pools
        output[cycle, "meanPanHet"] <- meanG(Hybrids) - mean(GP.EliteA$mu_HW, GP.EliteB$mu_HW)
        
        
        
      } #end cycle loop
      
      saveRDS(output, file.path(paste(scenario, "_",
                           "H2", 1, "_",
                           paste(paste0(colnames(scenParam), scenParam[1,]), collapse = "_"),
                           ".RDS",
                           sep = "")))
      
    #} #end REP loop
    
  #}# end START loop
  
#}# end H2 loop



