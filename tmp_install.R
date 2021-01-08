# from CRAN
install.packages(c("R.utils", "Rcpp", "here", "remotes", "devtools",
                   "data.table", "tidyverse", "origami", "hal9001",
                   "foreach", "doRNG", "future", "future.apply", "doFuture"))

# use remotes to install from GitHub
remotes::install_github(c("tlverse/sl3@master",
                          "nhejazi/haldensify@master",
                          "nhejazi/txshift@master"))
