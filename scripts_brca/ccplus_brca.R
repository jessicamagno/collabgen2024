#############################################################
###### BRCA Consensus cluster plus ##########################

#-- Loading packages
library(ConsensusClusterPlus)

#-- Loading BRCA LumB data
load("data_brca/BRCA_LumB_final_counts.RData")
load("data_brca/BRCA_LumB_final_clin.RData")

#-- Checking if counts and clinical data are in the same order
identical(colnames(gexp_LumB), rownames(clinical_BRCA_LumB)) #TRUE

d <- gexp_LumB
d <- sweep(d,1, apply(d,1,median,na.rm=T))

results <- ConsensusClusterPlus(d,reps=50,pItem=0.8,pFeature=1,clusterAlg="hc",
                               distance="pearson",plot="png", title = "brca_lumB_jul15")


#-- Generating cluster and item consensus
icl <- calcICL(results,plot="png")
icl[["clusterConsensus"]]

