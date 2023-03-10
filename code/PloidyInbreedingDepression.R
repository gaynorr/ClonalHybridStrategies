################################################################################
# This script demonstrates that with digenic dominance, autopolyploids have higher
# inbreeding depression values than their corresponding diploids generated by 
# genome reduction
################################################################################

library(AlphaSimR)

results <- data.frame(tet.HWE = NA,
                      tet.FullInbr = NA,
                      dip.HWE = NA,
                      dip.FullInbr = NA)

for(i in 1:10){
  # Create a founder population in LD with allele frequency 0.5
  founderPop = quickHaplo(100, #nInd
                          10,  #nChr
                          1000, #segSites
                          ploidy=4, #ploidy 
                          inbred=TRUE) #fully inbred founders
  
  
  # Add a trait with dominance
  SP = SimParam$
    new(founderPop)$
    addTraitAD(1000,meanDD=0.2,varDD=0.1)$
    setVarE(H2=1)
  
  
  # Create a tetraploid population, then reduce genome to get diploids
  tetraploid = newPop(founderPop)
  diploid = reduceGenome(tetraploid)
  
  GPT <- genParam(tetraploid)
  results[i, "tet.HWE"] <- GPT$mu_HW # tetraploid population mean at HWE
  results[i, "tet.FullInbr"] <- meanG(tetraploid) # tetraploid fully inbred population mean
  
  GPD <- genParam(diploid)
  results[i, "dip.HWE"] <- GPD$mu_HW # diploid population mean at HWE
  results[i, "dip.FullInbr"] <- meanG(diploid) # diploid fully inbred population mean
}

results <- round(results, digits = 2)
results$tet.InbrDep <- results$tet.FullInbr - results$tet.HWE
results$dip.InbrDep <- results$dip.FullInbr - results$dip.HWE

#write.csv(results, "C:/Users/MRLABROO/Dropbox (CIMMYT SeeD)/EiB Modules/Module 2/PMP/ResearchPapers/hybridBreeding/rawRevision/PloidyInbreedingDepression.csv")
