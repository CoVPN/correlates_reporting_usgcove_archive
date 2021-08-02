try.warning <- function(expr, silent=TRUE){

  # make it to error
  op <- options("warn")
  on.exit(options(op))
  options(warn = 2)

  # catch using the try() function

  try(expr, silent)
}

#
#try.error <- function(expr, silent=TRUE){
#
#  # make it to error
#  op <- options("warn")
#  on.exit(options(op))
#  options(warn = 1)
#
#  # catch using the try() function
#
#  try(expr, silent)
#}


try.error<-function(expr, silent=TRUE){

  # make it to error
  op <- options("warn")
  on.exit(options(op))
  options(warn = -1)

  # catch using the try() function

  try(expr, silent=TRUE)
}

resample <- function(x, ...) x[sample.int(length(x), ...)]

ee<-function(a,b){
abs(a-b)<1e-10
}

correct<-function(x,s)
{
  out<-NA
  if (!is.na(x)){
  if (x==0) out<-x
  if (x>0) out<-x-0.5*s
  if (x<0) out<-x+0.5*s
  }
  return(out)
}



##### distribution of linear combination of Y conditional on disease status

getAUC.linear<-function(mu.D,mu.Dbar,Sigma){
mu.c.D<-t(mu.D-mu.Dbar)%*%solve(Sigma)%*%mu.D
mu.c.Dbar<-t(mu.D-mu.Dbar)%*%solve(Sigma)%*%mu.Dbar

Sigma.c<-t(mu.D-mu.Dbar)%*%solve(Sigma)%*%(mu.D-mu.Dbar)

a=(mu.c.D-mu.c.Dbar)/sqrt(Sigma.c)
b=1

AUC=pnorm(a/sqrt(1+b^2))

return(AUC)

}

phid<-function(x)
{
 dnorm(x)*(-x)
}

pred.density<-function(x,predx){
 lapply(predx, function(ci) {
    out<-density(x,from=ci,to=ci,n=1)$y
    return(out)
    }
 )
}

xvl<-function(v,lambda)
{
 if (lambda==0) out<-qnorm(v)
 if (lambda!=0) out<-(exp(qnorm(v))^lambda-1)/lambda
 return(out)
}



xl<-function(x,lambda)
{
 if (lambda==0) out<-log(x)
 if (lambda!=0) out<-(x^lambda-1)/lambda
 return(out)
}

interROC<-function(ROC,tt){
 TPF<-ROC$sensitivity
 FPF<-1-ROC$specificity
 oo<-order(FPF)
 FPF<-FPF[oo]
 TPF<-TPF[oo]
 l1<-length(FPF)-N.G.E(tt,FPF)
 l2<-N.L.E(tt,FPF)
 FPF<-c(0,FPF)
 TPF<-c(0,TPF)
 l1<-l1+1
 l2<-l2+1
 ####### check whether tt falls on FPF point ########
 ### if yes, the TPF corresponds to the largest one with FPF=tt ###
 ### Reason: quantile function is defined to be left continuous, so goes to the smallest threshold  ####
 ### if no, interpolate ##########

 out<-rep(NA,length(tt))
 for (k in 1:length(tt)){
 out[k]<-ifelse(l1[k]!=l2[k],max(TPF[FPF==FPF[l2[k]]],na.rm=T),TPF[l1[k]]+(TPF[l1[k]+1]-TPF[l1[k]])/(FPF[l1[k]+1]-FPF[l1[k]])*(tt[k]-FPF[l1[k]]))
 }
 return(out)
}


#stpfun<-function(v,RV){
# out<-stepfun(v,c(0,RV),right=FALSE,f=1)
#}

stpfun<-function(v,RV){
 out<-stepfun(v,c(min(RV),RV),right=FALSE,f=1)
}



calU<-function(x){
 I11<-ifelse(x[1]>0,1,.5)
 I11<-ifelse(x[1]<0,0,I11)
 I22<-ifelse(x[2]>0,1,.5)
 I22<-ifelse(x[2]<0,0,I22)
 I12<-ifelse(x[1]+x[2]>0,1,.5)
 I12<-ifelse(x[1]+x[2]<0,0,I12)
 I13<-ifelse(x[1]+x[3]>0,1,.5)
 I13<-ifelse(x[1]+x[3]<0,0,I13)
 I34<-ifelse(x[3]+x[4]>0,1,.5)
 I34<-ifelse(x[3]+x[4]<0,0,I34)
 out<-c(I11,I22,I12,I13,I34)
 return(out)
}

calWU<-function(x,y){
 I11<-ifelse(x[1]<y[1],1,.5)
 I11<-ifelse(x[1]>y[1],0,I11)
 I22<-ifelse(x[2]<y[2],1,.5)
 I22<-ifelse(x[2]>y[2],0,I22)
 I12<-ifelse(x[1]<y[2],1,.5)
 I12<-ifelse(x[1]>y[2],0,I12)
 I13<-ifelse(x[1]<y[3],1,.5)
 I13<-ifelse(x[1]>y[3],0,I13)
 I34<-ifelse(x[3]<y[4],1,.5)
 I34<-ifelse(x[3]>y[4],0,I34)
 out<-c(I11,I22,I12,I13,I34)
 return(out)
}

calmoment<-function(data,w,nnx=n.Dbar,nny=n.D){
 ss<-data[w,]
 Y.Dbar<-ss[1:nnx,]
 Y.D<-ss[(nnx+1):(nnx+nny),]
 Y1.Dbar<-Y.Dbar[,1]
 Y2.Dbar<-Y.Dbar[,2]
 Y1.D<-Y.D[,1]
 Y2.D<-Y.D[,2]

 s1.p<-calpercentNM(Y1.Dbar,Y1.D)
 s2.p<-calpercentNM(Y2.Dbar,Y2.D)

 s1.e<-calpercentNP(Y1.Dbar,Y1.D)
 s2.e<-calpercentNP(Y2.Dbar,Y2.D)

 s.p<-s1.p-s2.p
 s.e<-s1.e-s2.e


 ####
 oo<-sample(n.D,replace=F)

 II.p1<-calU(s.p[oo[1:4]])
 WI.p1<-calWU(s1.p[oo[1:4]],s2.p[oo[1:4]])

 II.p2<-calU(s.p[oo[5:8]])
 WI.p2<-calWU(s1.p[oo[5:8]],s2.p[oo[5:8]])

 II.p3<-calU(s.p[oo[9:12]])
 WI.p3<-calWU(s1.p[oo[9:12]],s2.p[oo[9:12]])

 II.p4<-calU(s.p[oo[13:16]])
 WI.p4<-calWU(s1.p[oo[13:16]],s2.p[oo[13:16]])

 II.e1<-calU(s.e[oo[1:4]])
 WI.e1<-calWU(s1.e[oo[1:4]],s2.e[oo[1:4]])

 II.e2<-calU(s.e[oo[5:8]])
 WI.e2<-calWU(s1.e[oo[5:8]],s2.e[oo[5:8]])

 II.e3<-calU(s.e[oo[9:12]])
 WI.e3<-calWU(s1.e[oo[9:12]],s2.e[oo[9:12]])

 II.e4<-calU(s.e[oo[13:16]])
 WI.e4<-calWU(s1.e[oo[13:16]],s2.e[oo[13:16]])

 kk.e<-calstat(s1.e,s2.e)
 kk.p<-calstat(s1.p,s2.p)
 out<-c(kk.e,kk.p,II.e1,II.e2,II.e3,II.e4,WI.e1,WI.e2,WI.e3,WI.e4,II.p1,II.p2,II.p3,II.p4,WI.p1,WI.p2,WI.p3,WI.p4)
 return(out)
}



 calCI<-function(Est,SE,alpha=0.05,lowbound=NULL,highbound=NULL){
  level<-1-alpha/2
  low<-Est-qnorm(level)*SE
  if (!missing(lowbound)) low<-ifelse(low<lowbound,lowbound,low)
  high<-Est+qnorm(level)*SE
  if (!missing(highbound)) high<-ifelse(high>highbound,highbound,high)
  pvalue<-2*(1-pnorm(abs(Est/SE)))
  return(list(Est=Est,low=low,high=high,pvalue=pvalue))
 }

 calCI.logit<-function(Est,SE,alpha=0.05,mix=F){
  level<-1-alpha/2
  if (mix==F){
  if (is.na(SE)){
  pvalue<-low<-high<-NA
  }
  else if (!is.na(SE) & SE>0){
  SE<-SE/Est/(1-Est)
  Est<-log(Est/(1-Est))
  low<-Est-qnorm(level)*SE
  high<-Est+qnorm(level)*SE
  pvalue<-2*(1-pnorm(abs(Est/SE)))
  Est<-exp(Est)/(1+exp(Est))
  low<-exp(low)/(1+exp(low))
  high<-exp(high)/(1+exp(high))
  } else {
  pvalue<-NA
  Est<-Est
  low<-Est
  high<-Est
  }

  }
  if (mix==T){
  Est.o<-Est
  SE<-ifelse((ee(Est.o,0) | ee(Est.o,1)),SE,SE/Est/(1-Est))
  Est<-ifelse((ee(Est.o,0) | ee(Est.o,1)),Est,log(Est/(1-Est)))
  low<-ifelse(ee(Est.o,0), 0, Est-qnorm(level)*SE)
  high<-ifelse(ee(Est.o,1),1, Est+qnorm(level)*SE)
  pvalue<-2*(1-pnorm(abs(Est/SE)))
  Est<-ifelse((ee(Est.o,0) | ee(Est.o,1)),Est.o, exp(Est)/(1+exp(Est)))
  low<-ifelse((ee(Est.o,0) | ee(Est.o,1)),low,exp(low)/(1+exp(low)))
  high<-ifelse((ee(Est.o,0) | ee(Est.o,1)),high,exp(high)/(1+exp(high)))
  }
  return(list(Est=Est,low=low,high=high,pvalue=pvalue))
 }

 calCI.logit.new<-function(Est.logit,SE.logit,alpha=0.05,mix=F){
  level<-1-alpha/2
  low.logit<-Est.logit-qnorm(level)*SE.logit
  high.logit<-Est+qnorm(level)*SE.logit
  pvalue<-2*(1-pnorm(abs(Est/SE)))
  Est<-exp(Est.logit)/(1+exp(Est.logit))
  low<-exp(low.logit)/(1+exp(low.logit))
  high<-exp(high.logit)/(1+exp(high.logit))
  return(list(Est=Est,low=low,high=high,pvalue=pvalue))
 }
## count how many YY's are smaller or equal to yy

N.L.E <- function(yy, YY)  ## sum I(YY <= yy[i])
  {
    rank(c(yy+1e-8,YY))[1:length(yy)] - rank(yy)  ### add a small pertubation to avoid ties when calculating rank
  }

N.G.E <- function(yy, YY)  ## sum I(YY >= yy[i])
  {
    length(YY)-(rank(c(yy-1e-8,YY))[1:length(yy)] - rank(yy))
  }


#######


 hosmerlem <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0,
        1, 1/g)), include.lowest = T)
    obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    nk<-apply(expect,1,sum)
    chisq <- sum(((obs - expect)^2/(expect*(1-expect/nk)))[,2])
    P <- 1 - pchisq(chisq, g - 2)
    c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}

 hosmerlem <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0,
        1, 1/g)), include.lowest = T)
    obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    chisq <- sum((obs - expect)^2/expect)
    P <- 1 - pchisq(chisq, g - 2)
    c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}


 hosmerlem <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = unique(quantile(yhat, probs = seq(0,
        1, 1/g))), include.lowest = T)
    obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    chisq <- sum((obs - expect)^2/expect)
    gg<-length(unique(cutyhat))
    P <- 1 - pchisq(chisq, g - 2)
    c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}


hosmerlem.plot.old <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0,
        1, 1/g)), include.lowest = T,right=FALSE)
    obs <- xtabs(cbind(1-y,y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    ng<-apply(obs,1,sum)
    nk<-apply(expect,1,sum)
    RV.obs<-(obs/ng)[,2]
    RV.exp<-(expect/ng)[,2]
    return(list(obs=obs/nk,expect=expect/nk,v=seq(0,1,1/g),RV.obs=c(RV.obs[1],RV.obs),RV.exp=c(RV.exp[1],RV.exp)))
}
hosmerlem.plot <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0,
        1, 1/g)), include.lowest = T,right=FALSE)
    obs <- xtabs(cbind(1-y,y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    ng<-apply(obs,1,sum)
    nk<-apply(expect,1,sum)
    RV.obs<-(obs/ng)[,2]
    RV.exp<-(expect/ng)[,2]
    vv<-cumsum(as.numeric(table(cutyhat)/sum(table(cutyhat))))
    return(list(obs=obs/nk,expect=expect/nk,v=(c(0,vv[-g])+vv)/2,RV.obs=RV.obs,RV.exp=RV.exp))
}


hosmerlem.plot <-
function (y, yhat, g = 10)
{
    cutyhat <- cut(yhat, breaks = unique(quantile(yhat, probs = seq(0,
        1, 1/g))), include.lowest = T,right=FALSE)
    obs <- xtabs(cbind(1-y,y) ~ cutyhat)
    expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
    ng<-apply(obs,1,sum)
    nk<-apply(expect,1,sum)
    RV.obs<-(obs/ng)[,2]
    RV.exp<-(expect/ng)[,2]
    gg<-length(unique(cutyhat))
    vv<-cumsum(as.numeric(table(cutyhat)/sum(table(cutyhat))))
    return(list(obs=obs/nk,expect=expect/nk,v=(c(0,vv[-gg])+vv)/2,RV.obs=RV.obs,RV.exp=RV.exp))
}


bar.overlay<-function(data1,data2,label1,label2,label3,col1="lightblue",col2="mistyrose"){
hh1<-hist(data1,breaks=seq(0,1,.1),plot=F,prob=T)
 hh2<-hist(data2,breaks=seq(0,1,.1),plot=F,prob=T)
comb<-rbind(hh1$density,hh2$density)
row.names(comb)<-c(label1,label2)
colnames(comb)<-c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
barplot(comb,beside=T,legend=rownames(comb),col=c(col1,col2),ylab="density",xlab=label3)
#mtext(side = 2, at = mean(c(hh1$density,hh2$density))*1.5, line = 2,text = "density")
#mtext(side = 1, at = 10, line = 2,  text = "percentile values of CA-125")
}



bar.overlay4<-function(data1,data2,data3,data4,label1,label2,label3,label4,label5,col1="lightblue3",col2="lightcyan",
col3="lightpink3",col4="mistyrose"){
 hh1<-hist(data1,breaks=seq(0,100,10),plot=F,prob=T)
 hh2<-hist(data2,breaks=seq(0,100,10),plot=F,prob=T)
 hh3<-hist(data3,breaks=seq(0,100,10),plot=F,prob=T)
 hh4<-hist(data4,breaks=seq(0,100,10),plot=F,prob=T)
comb<-rbind(hh1$density,hh2$density,hh3$density,hh4$density)
row.names(comb)<-c(label1,label2,label3,label4)
colnames(comb)<-c("0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100")
barplot(comb,beside=T,legend=rownames(comb),col=c(col1,col2,col3,col4),
ylab="density",xlab=label5)
#mtext(side = 2, at = mean(c(hh1$density,hh2$density))*1.5, line = 2,text = "density")
#mtext(side = 1, at = 10, line = 2,  text = "percentile values of CA-125")
}


boundary.p <- function(densk) {
        for(i in 1:length(densk$x)){
            if(densk$x[i]<0) {
                i1 <- which.max(densk$x>abs(densk$x[i]))
                densk$y[i1] <- densk$y[i] + densk$y[i1]
                densk$y[i]<-0
            }
            if(densk$x[i]>100) {
                i2 <- which.min(densk$x<(200-densk$x[i]))
                densk$y[i2] <- densk$y[i] + densk$y[i2]
                densk$y[i]<-0
            }
        }
        dat1 <- data.frame(densk$x, densk$y)
        s1 <- subset(dat1, (dat1$densk.x>=0)&(dat1$densk.x<=100))
        return(s1)
    }

#fit00<-lrm(highgrade~lpsa, data=riskdata, x=TRUE, y=TRUE)
#fit11<-lrm(highgrade~lpsa+age+dre+priorbiop,data=riskdata,x=TRUE, y=TRUE)
#fit22<-lrm(highgrade~simY,data=riskdata,x=TRUE, y=TRUE)



#residuals.lrm(fit00,type='gof')
#residuals.lrm(fit11,type='gof')
#residuals.lrm(fit22,type='gof')
