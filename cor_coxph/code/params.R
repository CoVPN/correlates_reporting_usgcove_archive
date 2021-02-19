study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well

# 1e3 bootstraps take 8 min with 30 CPUS. The results are saved in several .Rdata files.
# permutation also takes time
# DB: maybe parallel::detectCores() better than 30?
numCores=if(Sys.info()['sysname']=="Windows") 1 else 30 # number of cores available on the machine
# B=1000 # number of bootstrap replicates
B=25 # number of bootstrap replicates

# DB: Should seeds be placed in params.R?
# seeds=1:1e4
seeds=1:25

#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")

trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")
