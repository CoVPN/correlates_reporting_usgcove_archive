numCores <- unname(ifelse(Sys.info()["sysname"] == "Windows",
                          1, future::availableCores()))

#if (length(assays) %in% c(3,4)) {
#  .mfrow <- c(2, 2)
#} else if (length(assays) == 5) {
#  .mfrow <- c(3, 2)
#} else if (length(assays) == 2) {
#  .mfrow <- c(1, 2)
#} else {
#  stop("Please re-define variable .mfrows")
#}
.mfrow <- c(1, 1)
