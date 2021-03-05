study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for D57 CoR
B=1e3 # number of bootstrap replicates 1e3
numPerm=1e4 # number permutation replicates 1e4
numCores=if(Sys.info()['sysname']=="Windows") 1 else parallel::detectCores()

.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")
