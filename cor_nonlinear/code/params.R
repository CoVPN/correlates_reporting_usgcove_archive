study.name="mock" # study.name is used in figure/table file names and printed in tables/figures as well
.mfrow=if(length(assays)==4) c(2,2) else if(length(assays)==5) c(3,2) else stop("pls redefine .mfrows")
