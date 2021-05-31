if (length(assays) %in% c(3,4)) {
  .mfrow <- c(2, 2)
} else if (length(assays) == 5) {
  .mfrow <- c(3, 2)
} else {
  stop("Please re-define variable .mfrows")
}
