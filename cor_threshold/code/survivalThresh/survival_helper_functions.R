
bound <- function(x, bounds) {
  lower <- bounds[[1]]
  if (length(bounds) > 1) {
    upper <- bounds[[2]]
  } else {
    upper <- 1 - lower
  }
  pmin(pmax(x, lower), upper)
}

#' @import simcausal
simulate <- function(n, censoring = TRUE) {

  do <- function(x) {
    x <- x
    x <- floor(x*5)
    x <- pmin(x,15)
    x <- pmax(x,1)

  }
  round <- Vectorize(round)

  D <- DAG.empty()
  D <- D +
    node("W1", distr = "runif", min = -1, max = 1 ) +
    node("W2", distr = "runif", min = -1, max = 1 ) +

    node("A", distr = "rbeta", shape1 = 3 +1.5*(W1 + W2), shape2 = 3 + 1.5*(W1 + W2) ) +
    node("J", distr = "rbeta", shape1 = 1.5 + t/15 + (W1 + W2 + A)/2, shape2 = 2, t = 1:15) +
    #node("J", distr = "rbeta", shape1 = 2, shape2 = 2, t = 1:15) +

    node("Nt", distr = "rbinom", size = 1, prob = (0.07-A/15+(t-1)/12)*plogis( t/15 -2*A+ (W1 + W2)/2), t = 1:15, EFU = TRUE)+
    #node("Ct", distr = "rbinom", size = 1, prob = censoring*(0.02)*plogis(t/15), t = 1:15, EFU = TRUE)
    node("Ct", distr = "rbinom", size = 1, prob = censoring*(0.03+(t-1)/40)*plogis(t/15 - A + (W1 + W2)/2), t = 1:15, EFU = TRUE)



  setD <- set.DAG(D)
  suppressMessages(dat <- sim(setD, n = n))


  dat_C <- dat[,grep("Ct_[0-9]+", colnames(dat), value = T)]
  dat_N <- dat[,grep("Nt_[0-9]+", colnames(dat), value = T)]
  dat_J <- dat[,grep("J_[0-9]+", colnames(dat), value = T)]
  .C <- apply(dat_C,1, function(Ct) {
    t <- which(!is.na(Ct) & Ct==1)
    if(length(t) == 0){
      t <- 14
    } else {
      t <- t[1]
    }
    return(t)
  })

  .T <- apply(dat_N,1, function(Nt) {
    t <- which(!is.na(Nt) & Nt==1)
    if(length(t) == 0){
      t <- 15
    } else {
      t <- t[1]
    }
    return(t)
  })

  .J <- apply(dat_J,1, function(Jt) {
    J <- Jt[max(which(!is.na(Jt)))]
    return(J)
  })

  data <- data.table(T = .T, C = .C, J = .J, W1 = dat$W1, W2 = dat$W2, A = dat$A, id= dat$ID)
  data$Ttilde <- pmin(data$T, data$C)
  data$Delta <- as.numeric(data$T <= data$C)
  return(data)

}

shape_long <- function(data, node_list, times = 1:max(data$Ttilde)) {
  times <- 1:max(data$Ttilde)
  t_mat <- matrix(rep(times, nrow(data)), ncol = length(times), byrow = T)
  Tt <- 1*(t_mat >= data$Ttilde)
  Nt <- data$Delta*Tt
  at_risk  <- 1*(t_mat <= data$Ttilde)
  Ct  <- (1-data$Delta)*Tt
  #at_risk_Ct <- at_risk_Nt*(1-Nt)
  data_long <- data.table(t = as.vector(t_mat), Nt = as.vector(Nt), Ct = as.vector(Ct), at_risk = as.vector(at_risk),id = rep(1:nrow(data), length(times)), Event = as.vector(1*(t_mat == data$Ttilde)))
  set(data_long, , node_list$W, data.table(apply(data[,node_list$W, with = F],2, rep, length(times))))
  set(data_long, , "A", rep(data[[node_list$A]],length(times)))
  set(data_long, , "J", rep(data[[node_list$J]],length(times)))
  set(data_long, , node_list$weights, rep(data[[node_list$weights]],length(times)))
  return(data_long)
}




hazard_to_survival <- function(v, nt, left = F) {
  s <- t(apply(matrix(1-v, ncol = nt),1,cumprod))
  if(left) {
    s <- cbind(rep(1, nrow(s)),s[,-ncol(s)])
  }
  return(as.vector(s))
}

compute_survival_functions <- function(likelihoods, nt) {
  Ft_list <- list()
  surv_N <- hazard_to_survival(likelihoods$N, nt = nt, left = T)
  Ft_p1 <- matrix(surv_N, ncol = nt)*matrix(likelihoods$N, ncol = nt)
  for(i in 1:ncol(likelihoods$J)) {
    Ft <- Ft_p1*matrix(likelihoods$J[[i]], ncol = nt)
    Ft <- t(apply(Ft,1, cumsum))
    Ft_list[[i]] <- as.vector(Ft)

  }
  Ft <- as.matrix(do.call(cbind, Ft_list))
  rm(Ft_list)
  return(list(St = surv_N, Ft = Ft))
}

