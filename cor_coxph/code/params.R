study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS. The results are saved in several .Rdata files.
B=1e3 # number of bootstrap replicates
numPerm=1e4 # number permutation replicates
numCores=if(Sys.info()['sysname']=="Windows") 1 else parallel::detectCores()

#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")

trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")
