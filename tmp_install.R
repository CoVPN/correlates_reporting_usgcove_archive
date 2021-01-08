# from CRAN
install.packages(c("R.utils", "Rcpp", "here", "remotes", "devtools",
                   "data.table", "tidyverse", "origami", "hal9001",
                   "foreach", "doRNG", "future", "future.apply", "doFuture"))

# use remotes to install from GitHub
remotes::install_github(c("tlverse/sl3@master",
                          "nhejazi/haldensify@master",
                          "nhejazi/txshift@master"))

if (!("osDesign" %in% installed.packages())) {
  pkg_file <- "osDesign_1.7.tar.gz"
  url <- paste0("http://cran.r-project.org/src/contrib/Archive/osDesign/",
                pkg_file)
  download.file(url = url, destfile = pkg_file)
  install.packages(pkgs = pkg_file, type = "source", repos = NULL)
  unlink(pkg_file)
}

## install local data package
remotes::install_github("youyifong/CovidCorrSAP/R_packages/COVIDcorr")
#remotes::install_local(here("..", "R_packages", "COVIDcorr"), force = TRUE)
library(COVIDcorr); stopifnot(packageVersion("COVIDcorr") >= "2020.10.14")
