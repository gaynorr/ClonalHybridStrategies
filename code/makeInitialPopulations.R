library(AlphaSimR)
#Rcpp::sourceCpp("GPCP_polyploids.cpp")




meanDD <- c(0, 0.5, 1, 1.5, 10) #mean dominance degree
varDD <- c(0, 0.2, 1, 10) #variance of dominance degrees
#H2 <- c(0, 0.25, 0.5, 0.75, 1) #broad-sense heritability- doesn't need a starter, can be set in scenario
ploidy <- c(2, 4, 6)
nQTL <- c(100, 1000, 5000) #number of QTL affecting trait
nSNP <- 100 #c(100, 1000, 5000) #number of SNPs genotyped
nChr <- 10
nReps <- 10

# Create one population per rep per architecture

starters <- expand.grid(list(meanDD, varDD, ploidy, nQTL, nSNP, nChr)) #cross combinations
colnames(starters) <- c("meanDD", "varDD", "ploidy", "nQTL", "nSNP", "nChr")

for(i in 1:nrow(starters)){
  for(j in 1:nReps){
  print(i)
  
  founderPop <- runMacs(nInd = 80,
                        nChr = starters[i, "nChr"],
                        segSites = starters[i, "nQTL"] + starters[i, "nSNP"],
                        ploidy = starters[i , "ploidy"])
  
  SP = SimParam$
    new(founderPop)$
    addTraitAD(nQtlPerChr = starters[i, "nQTL"],
               mean = 0,
               var = 1,
               meanDD = starters[i, "meanDD"],
               varDD = starters[i, "varDD"])
  SP$addSnpChip(nSnpPerChr = starters[i, "nSNP"])
  
  scenParam <- starters[i, ]
  scenParam$rep <- j
  
  outname <- paste(paste(colnames(starters), starters[i,], sep = "", collapse = "_"),  
                   paste("_rep", j, sep = ""), sep = "", collapse = "")
  save(founderPop, scenParam, SP, file = file.path(paste(outname, ".RData", sep = "")))
  }
}
