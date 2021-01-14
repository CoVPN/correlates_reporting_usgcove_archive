# Function to truncate at LLOQ and ULOQ
setLOQ <- function(x, lloq=-Inf, uloq=Inf){
  lloq <- ifelse(is.infinite(lloq), lloq, log10(lloq))
  uloq <- log10(uloq)
  x <- ifelse(x < lloq, lloq-log10(2), x)
  x <- ifelse(x > uloq, uloq, x)
}

# Function to Generate Responder calls
setResponder <- function(data, bl, post, folds=c(2, 4), lloq, responderFR=4){
  data[, paste0(post, "FR")] <- 10^(data[, post]-data[, bl])
  foldsInd <- sapply(folds, function(x) as.numeric(data[, paste0(post, "FR")] >= x)) %>% data.frame()
  names(foldsInd) <- paste0(post, "FR", folds)
  data <- cbind(data, foldsInd)
  data[, paste0(post, "Resp")] <- ifelse((data[, bl] < log10(lloq) & data[, post]>log10(lloq)) |
                                         (data[, bl] >=log10(lloq) & data[, paste0(post, "FR", responderFR)]==1), 
                                         1, 
                                         0)
  return(select_at(data, c(paste0(post, "Resp"), paste0(post, "FR", folds))))
}

# Function Generate delta between T1 and T2
setDelta <- function(data, endpoint, t1, t2, prefix="", suffix=""){
  data[, paste0(prefix, endpoint, "delta", suffix)] <- data[paste0(t2, endpoint)] - data[paste0(t1, endpoint)]
  return(select_at(data, paste0(prefix, endpoint, "delta", suffix)))
}

# Function to Generate summaries for numeric demographic variables or categorical variables
demoSumm <- function(data, numVar, catVar){
  retNum <- map_dfr(data[, numVar], .id="covariate", function(Subgroup){
    iCat <- summary(Subgroup)
    if (class(Subgroup)=="numeric") iCat <- round(iCat, 1)
    iCat
    })
  retCat <- map_dfr(data[, catVar], .id="covariate", function(Subgroup){
    xlen <- length(Subgroup)
    iCat <- table(Subgroup) %>% data.frame() %>% mutate(Prop=Freq/xlen, Pct=sprintf("%.1f%%", 100*Prop))
    iCat
  })
  return(list(numeric=retNum, category=retCat))
}
