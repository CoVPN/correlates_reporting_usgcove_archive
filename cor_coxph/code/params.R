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


# need this function b/c svycoxh may error due to singularity if, e.g. all cases have the same marker value
run.svycoxph=function(f, design) {
    fit=try(svycoxph(f, design=design), silent=T)
    if (class(fit)[1]=="try-error") NA else fit
}


# 26Oct2020      Erika Rudnicki
theforestplot <- function(cohort=NA,group,nEvents=NA,totFU=NA,rate=NA,point.estimates,lower.bounds,upper.bounds,p.values,table.labels,zero.line=1.0,dashed.line=NA,
    x.ticks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
    decimal.places = 1,fontsize = 1,width.pdf=7,height.pdf=7,graphwidth="auto",...){
  
  plotdata <- structure(
    list(
      mean  = c(NA, point.estimates),
      lower = c(NA, lower.bounds),
      upper = c(NA, upper.bounds)
    ),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA,-(length(point.estimates) + 1)),
    class = "data.frame"
  )
  
  # remove.redundancy <- function(x){ifelse(!is.na(dplyr::lag(x)) & x==dplyr::lag(x), NA, x)}
  show.decimals <- function(x){format(round(x, decimal.places), nsmall=decimal.places)}
  
  group_edit <- sapply(group, function(x){
    if(grepl(">=", x)){
      ssplit <- strsplit(x, split = ">=")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" >= "', ssplit[2], '")')))
    }else if(grepl("<", x)){
      ssplit <- strsplit(x, split = "<")[[1]]
      eval(parse(text = paste0('expression("', ssplit[1], '" < "', ssplit[2], '")')))
    }else{
      x
    }
  }, simplify = FALSE)
  group <- group_edit
    
  if(all(is.na(p.values))){
    tabletext <- list(
      # c(table.labels[1], remove.redundancy(as.character(cohort))),
      c(table.labels[1], group),
      c(table.labels[3], nEvents),
      # c(table.labels[4], totFU),
      # c(table.labels[5], rate),
      c(paste0(table.labels[2]), 
        paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
      c(" ", rep(NA, length(point.estimates)))
    )
  } else{
      tabletext <- list(
        # c(table.labels[1], remove.redundancy(as.character(cohort))),
        c(table.labels[1], group),
        c(table.labels[3], nEvents),
        # c(table.labels[4], totFU),
        # c(table.labels[5], rate),
        c(paste0(table.labels[2]), 
          paste0(sapply(point.estimates, show.decimals), " (", sapply(lower.bounds, show.decimals), ", ", sapply(upper.bounds, show.decimals), ")")),
        c("P-value", p.values)
      )}    
    
  replaceNA <- function(x){ifelse(grepl("NA", x), NA, x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceNA)
  
  replaceDash <- function(x){gsub("-", "\u2013", x)}
  tabletext[[3]] <- sapply(tabletext[[3]], replaceDash)
  
  if(!is.na(dashed.line)){grid.line <- structure(dashed.line, gp = gpar(lty = 2, col = "red", lwd=0.5))} else{grid.line <- FALSE}
  
  forestplot(tabletext,plotdata,is.summary = FALSE,col = fpColors(box = "darkblue",line = "darkblue",summary = "royalblue",zero="black"),    
    graph.pos = 3,graphwidth = graphwidth,hrzl_lines = list("2" = gpar(lty=1)),zero = zero.line,lwd.zero = 0.5,lwd.ci = 0.5,lwd.xaxis = 0.5,xticks = x.ticks,boxsize = 0.1,grid=grid.line,txt_gp = fpTxtGp(
      ticks = gpar(fontfamily = "", cex = fontsize * 0.8),
      label = gpar(fontfamily = "", cex = fontsize * 0.9),
      summary = gpar(cex = fontsize)
    ),
    colgap = unit(2, "mm"),align = c("l", "l", "l"),mar = unit(c(4,1,9,1), "mm"), #bltr
    clip = c(min(x.ticks), max(x.ticks)), ...
  )
}


# for some subpopulations, some strata may have empty ph2 
# this function removes those strata
get.dat.with.no.empty=function(dat.tmp) {
    tab=with(dat.tmp, table(Wstratum, ph2))
    subset(dat.tmp, !Wstratum %in% as.integer(rownames(tab)[which(tab[,"TRUE"]==0)]))
}
