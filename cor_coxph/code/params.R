study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well

# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS. The results are saved in several .Rdata files.
B=1e3 # number of bootstrap replicates 1e3
numPerm=1e4 # number permutation replicates 1e4
numCores=if(Sys.info()['sysname']=="Windows") 1 else parallel::detectCores()

#assays=c("bindSpike","bindRBD","pseudoneutid50","liveneutmn50","pseudoneutid80")
assays=c("bindSpike","bindRBD","pseudoneutid50","pseudoneutid80")
.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")

trt.labels=c("Placebo","Vaccine")
bstatus.labels=c("Baseline Neg","Pos")

Bstratum.labels <- c(
  "Age >= 65",
  "Age < 65, At risk",
  "Age < 65, Not at risk"
)
names(Bstratum.labels) <- Bstratum.labels
llods=c(bindSpike=20, bindRBD=20, pseudoneutid50=10, pseudoneutid80=10, liveneutmn50=62)
