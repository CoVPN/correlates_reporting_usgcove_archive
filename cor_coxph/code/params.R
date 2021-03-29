# B=1e3 and numPerm=1e4 take 10 min to run with 30 CPUS for D57 CoR
B <- 5 # number of bootstrap replicates 1e3
numPerm <- 5 # number permutation replicates 1e4
numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows",
                          1, future::availableCores()))

if (length(assays) == 4) {
  .mfrow <- c(2, 2)
} else if (length(assays) == 5) {
  .mfrow <- c(3, 2)
} else {
  stop("Please re-define variable .mfrows")
}
