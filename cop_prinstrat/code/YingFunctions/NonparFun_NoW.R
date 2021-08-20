### here we assume W has no effect on disease risk after conditional on S  #####

EM.cc.CPV.Probit.Small.Nonpar<-function(Z,Dens,S1,W,Y,delta,IC,beta,smooth=FALSE,lambda=1,arg=1){

#### smooth: indiator of whether error for location scale family is function of an extra covaraite
#### lambda: lambda for kernel function used for smoothing
#### Ucov: covariate for loc parameter
#### Wcov: covariate for log(scale) parameter

if (arg==1){

ind<-rep(1:sum(delta==0 & IC==0),dim(Dens)[2])  ### for all these subjects in placebo arm and not in IC

YY<-rep(Y[delta==0 & IC==0],dim(Dens)[2])

#S1.res<-rep(S1.out$res,rep(sum(delta==0 & IC==0),sum(delta==1)))


#resWeights<-rep(S1.out$weights,rep(sum(delta==0 & IC==0),sum(delta==1))) ### weights since residuals are not from iid
WW<-rep(W[delta==0 & IC==0],dim(Dens)[2])
ZZ<-rep(Z[delta==0 & IC==0],dim(Dens)[2])



#Loc1<-Ucov[delta==0 & IC==0,]%*%S1.out$eta
#Scale1<-as.matrix(Wcov[delta==0 & IC==0,])%*%S1.out$delta
#
#LLoc1<-rep(Loc1,sum(delta==1))
#SScale1<-rep(Scale1,sum(delta==1))
#SS1<-S1.res*exp(SScale1)+LLoc1

SS1<-as.numeric(rep(colnames(Dens),rep(sum(delta==0 & IC==0), dim(Dens)[2])))





 oo<-match(W[delta==0 & IC==0],rownames(Dens))
 resWeights<-as.numeric(as.matrix(Dens[oo,]))

   # if (smooth==TRUE) {
#        WW.val<-rep(W[delta==1],rep(sum(delta==0 & IC==0),sum(delta==1)))
#        smoothWeights<-exp(-0.5*((WW-WW.val)/lambda)^2)
#        denom.smoothWeights<-tapply(smoothWeights,ind,sum)
#        denom.smoothWeights<-rep(denom.smoothWeights,sum(delta==1))
#        smoothWeights<-smoothWeights/denom.smoothWeights
#
#    }
#
}


beta.old<-beta

#################### arg=2 is not updated yet compared to EM.cc ############################
if (arg==2){

n.C<-sum(Z==1 | IC==1) ### number of complete data
n.R<-sum(IC==0 & Z==0) ### number of reduced data

ind<-rep(1:sum(IC==0 & Z==0),sum(Z==1 | IC==1))  ### for all these subjects in placebo arm and not in IC

ind.Z<-rep(1:n.C,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1))) ### indicator for residual

YY<-rep(Y[IC==0 & Z==0],sum(Z==1 | IC==1))

oo<-IC==1
S1.res.IC<-(S1[oo]-S1.out$eta[1]-S1.out$eta[2]*W[oo])/exp(S1.out$delta)


S1.res<-rep(c(S1.out$res,S1.res.IC),rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))

WW<-rep(W[IC==0 & Z==0],sum(Z==1 | IC==1))
SS1<-S1.res*(exp(S1.out$delta))+S1.out$eta[1]+S1.out$eta[2]*WW


g.old<-rep(1/n.C,n.R)

}
iter<-0
repeat{

### probability of delta=1 conditional on S1,W

### if random sample S from Z=1
#P.Z1<-mean(delta[Z==1 & Y==1])*mean(Y[Z==1])+mean(delta[Z==1 & Y==0])*mean(1-Y[Z==1])

### if random sample S from cases and controls in vaccine arm
#P.Z1<-mean(delta[Z==1 & Y==1])*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+mean(delta[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))


P.Z1<-mean(delta[Z==1])

P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
ff<-ff*resWeights
ff<-ff/PP
if (smooth==TRUE){
    ff<-ff*smoothWeights
}


if (arg==1) {
 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)

 denom<-rep(denom,sum(delta==1))
 weights<-ff/denom
 } else weights<-ff

}  else {

Ljj.old<-tapply(ff*g.old,ind,sum)
Ljj.old<-rep(Ljj.old,sum(Z==1 | IC==1))
g.new<-tapply(ff*g.old/Ljj.old,ind.Z,sum)
g.new<-rep(g.new,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))
g.new<-(1+g.new)/(n.C+n.R)

if (length(ff)>0) {

denom<-tapply(ff*g.new,ind,sum)
denom<-rep(denom,sum(Z==1 | IC==1))
weights<-ff*g.new/denom
} else weights<-ff

}
#YY<-WW<-SS1<-SS2<-weights<-NULL

dataout<-data.frame(Y=c(Y[delta==1 | IC==1],YY),Z=c(Z[delta==1 | IC==1],ZZ),W=c(W[delta==1 | IC==1],WW),
S1=c(S1[delta==1 | IC==1],SS1),weights=c(rep(1,sum(delta==1 | IC==1)),weights))

#date()
#fit<-glm(Y~Z+S1+I(Z*S1)+W+I(Z*W),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
#fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>70 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>70 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}
############33


EM.cc.CPV.Probit.Small.Nonpar2<-function(Z,Sout,S1,W,Y,Wu,delta,IC,beta,smooth=FALSE,lambda=1,arg=1){

#### smooth: indiator of whether error for location scale family is function of an extra covaraite
#### lambda: lambda for kernel function used for smoothing
#### Ucov: covariate for loc parameter
#### Wcov: covariate for log(scale) parameter
Sval<-Sout$Sval
lval<-Sout$lval
Sweights<-Sout$Sweights
if (arg==1){


oo<-match(W[delta==0 & IC==0],Wu)
#ind<-rep(1:sum(delta==0 & IC==0),dim(Dens)[2])  ### for all these subjects in placebo arm and not in IC
ind<-rep(1:sum(delta==0 & IC==0), unlist(lval[oo]))

#YY<-rep(Y[delta==0 & IC==0],dim(Dens)[2])
YY<-rep(Y[delta==0 & IC==0],unlist(lval[oo]))

#S1.res<-rep(S1.out$res,rep(sum(delta==0 & IC==0),sum(delta==1)))


#resWeights<-rep(S1.out$weights,rep(sum(delta==0 & IC==0),sum(delta==1))) ### weights since residuals are not from iid
#WW<-rep(W[delta==0 & IC==0],dim(Dens)[2])
WW<-rep(W[delta==0 & IC==0],unlist(lval[oo]))


#ZZ<-rep(Z[delta==0 & IC==0],dim(Dens)[2])
ZZ<-rep(Z[delta==0 & IC==0],unlist(lval[oo]))


#Loc1<-Ucov[delta==0 & IC==0,]%*%S1.out$eta
#Scale1<-as.matrix(Wcov[delta==0 & IC==0,])%*%S1.out$delta
#
#LLoc1<-rep(Loc1,sum(delta==1))
#SScale1<-rep(Scale1,sum(delta==1))
#SS1<-S1.res*exp(SScale1)+LLoc1

#SS1<-as.numeric(rep(colnames(Dens),rep(sum(delta==0 & IC==0), dim(Dens)[2])))


SS1<-unlist(Sval[oo])

resWeights<-unlist(Sweights[oo])

 #oo<-match(W[delta==0 & IC==0],rownames(Dens))
 #resWeight<-Dens[oo,]
 #resWeight<-as.numeric(t(Dens[oo,]))

   # if (smooth==TRUE) {
#        WW.val<-rep(W[delta==1],rep(sum(delta==0 & IC==0),sum(delta==1)))
#        smoothWeights<-exp(-0.5*((WW-WW.val)/lambda)^2)
#        denom.smoothWeights<-tapply(smoothWeights,ind,sum)
#        denom.smoothWeights<-rep(denom.smoothWeights,sum(delta==1))
#        smoothWeights<-smoothWeights/denom.smoothWeights
#
#    }
#
}


beta.old<-beta

#################### arg=2 is not updated yet compared to EM.cc ############################
if (arg==2){

n.C<-sum(Z==1 | IC==1) ### number of complete data
n.R<-sum(IC==0 & Z==0) ### number of reduced data

ind<-rep(1:sum(IC==0 & Z==0),sum(Z==1 | IC==1))  ### for all these subjects in placebo arm and not in IC

ind.Z<-rep(1:n.C,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1))) ### indicator for residual

YY<-rep(Y[IC==0 & Z==0],sum(Z==1 | IC==1))

oo<-IC==1
S1.res.IC<-(S1[oo]-S1.out$eta[1]-S1.out$eta[2]*W[oo])/exp(S1.out$delta)


S1.res<-rep(c(S1.out$res,S1.res.IC),rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))

WW<-rep(W[IC==0 & Z==0],sum(Z==1 | IC==1))
SS1<-S1.res*(exp(S1.out$delta))+S1.out$eta[1]+S1.out$eta[2]*WW


g.old<-rep(1/n.C,n.R)

}
iter<-0
repeat{

### probability of delta=1 conditional on S1,W

### if random sample S from Z=1
#P.Z1<-mean(delta[Z==1 & Y==1])*mean(Y[Z==1])+mean(delta[Z==1 & Y==0])*mean(1-Y[Z==1])

### if random sample S from cases and controls in vaccine arm
#P.Z1<-mean(delta[Z==1 & Y==1])*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+mean(delta[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))


P.Z1<-mean(delta[Z==1])

P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
ff<-ff*resWeights
#ff<-ff/PP
if (smooth==TRUE){
    ff<-ff*smoothWeights
}


if (arg==1) {
 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)

 #denom<-rep(denom,sum(delta==1))
 denom<-denom[ind]
 weights<-as.numeric(ff)/denom
 } else weights<-ff

}  else {

Ljj.old<-tapply(ff*g.old,ind,sum)
Ljj.old<-rep(Ljj.old,sum(Z==1 | IC==1))
g.new<-tapply(ff*g.old/Ljj.old,ind.Z,sum)
g.new<-rep(g.new,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))
g.new<-(1+g.new)/(n.C+n.R)

if (length(ff)>0) {

denom<-tapply(ff*g.new,ind,sum)
denom<-rep(denom,sum(Z==1 | IC==1))
weights<-ff*g.new/denom
} else weights<-ff

}
#YY<-WW<-SS1<-SS2<-weights<-NULL

dataout<-data.frame(Y=c(Y[delta==1 | IC==1],YY),Z=c(Z[delta==1 | IC==1],ZZ),W=c(W[delta==1 | IC==1],WW),
S1=c(S1[delta==1 | IC==1],SS1),weights=c(rep(1,sum(delta==1 | IC==1)),weights))

#date()
#fit<-glm(Y~Z+S1+I(Z*S1)+W+I(Z*W),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
#fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>200 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>200 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}




##### allows case-control sampling of S among vaccinees also

Ssub.cc<-function(datain,ratio.v,ratio.p){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


if (ratio.p==-1){
    p.IC.1=1
} else {
p.IC.1<-ratio.p*sum(Y[Z==0]==1)/sum(Y[Z==0]==0)
}

p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)

IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==1) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


if (ratio.v==-1){
 p.CC.1=1
} else {
p.CC.1<-ratio.v*sum(Y[Z==1]==1)/sum(Y[Z==1]==0)
}
p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

#### do an extra sampling just to keep the seed same #######
if (p.CC==1) rbinom(N,1,0.1)

weights<-ifelse(Z==0,0,1)
weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
IC.vac<-IC.vac | (Z==1 & Y==1)


################# select subcohort for IC ######

datain$IC<-IC
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}


##### simulate data ###############

Ssub.old<-function(datain,p.IC){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)

################# select subcohort for IC ######

#IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)
#IC<-(Z==0)*rbinom(N,1,p.IC)
if (p.IC==1) IC<-(Z==0) & (Y==0)
else if (p.IC==0) IC<-rep(0,length(Z))
else IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

datain$IC<-IC
datain$IC.vac<-ifelse(Z==1,1,0)
datain$Y<-Y

weights<-ifelse(Z==0,0,1)
datain$weights<-weights
return(dataout=datain)
}

##### allows case-control sampling of S among vaccinees also

Ssub.cc<-function(datain,ratio.v,ratio.p){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


if (ratio.p==-1){
    p.IC.1=1
} else {
p.IC.1<-ratio.p*sum(Y[Z==0]==1)/sum(Y[Z==0]==0)
}

p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)

IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==1) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


if (ratio.v==-1){
 p.CC.1=1
} else {
p.CC.1<-ratio.v*sum(Y[Z==1]==1)/sum(Y[Z==1]==0)
}
p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

#### do an extra sampling just to keep the seed same #######
if (p.CC==1) rbinom(N,1,0.1)

weights<-ifelse(Z==0,0,1)
weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
IC.vac<-IC.vac | (Z==1 & Y==1)


################# select subcohort for IC ######

datain$IC<-IC
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}


##### simulate data ###############

Ssub.old<-function(datain,p.IC){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)

################# select subcohort for IC ######

#IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)
#IC<-(Z==0)*rbinom(N,1,p.IC)
if (p.IC==1) IC<-(Z==0) & (Y==0)
else if (p.IC==0) IC<-rep(0,length(Z))
else IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

datain$IC<-IC
datain$IC.vac<-ifelse(Z==1,1,0)
datain$Y<-Y

weights<-ifelse(Z==0,0,1)
datain$weights<-weights
return(dataout=datain)
}

##### allows random sampling of S among vaccinees also

Ssub<-function(datain,ratio.v,ratio.p){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)

p.IC.1<-ratio.p*sum(Y[Z==0]==1)/sum(Y[Z==0]==0)
p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)
IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==1) rbinom(N,1,0.1)

### simple random sampling of vaccinees ############
#p.CC.1<-ratio.v*sum(Y[Z==1]==1)/sum(Y[Z==1]==0)
#p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-(Z==1)*rbinom(N,1,ratio.v)
#IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

#### do an extra sampling just to keep the seed same #######
if (ratio.v==1) rbinom(N,1,0.1)


weights<-ifelse(Z==0,0,sum(Z==1)/sum(IC.vac))
#weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
#IC.vac<-IC.vac | (Z==1 & Y==1)


################# select subcohort for IC ######

datain$IC<-IC
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}



######
#
#R.probit<-function(Z,S1,W,beta){
######## calculate the logit risk ############
#   if (length(S1)>0) {
#   #Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z,W,W*Z)
#   #Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z,W)
#   Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z)
#   RR<-Xmatrix%*%beta
#   #Risk<-expit(RR)
#   Risk<-pnorm(RR)
#   } else Risk<-NULL
#
#   return(Risk)
#}

R.probit<-function(Z,S1,W,beta){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0) {
       Risk <- pnorm(cbind(1,Z,S1,S1*Z) %*% as.matrix(beta))
       gc()
   }
   return(Risk)
}

R.probit.L<-function(Z,S1,W,beta){
####### calculate the logit risk ############
   if (length(S1)>0) {
   #Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z,W,W*Z)
   #Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z,W)
   Xmatrix<-cbind(rep(1,length(S1)),Z,S1,S1*Z,W)
   RR<-Xmatrix%*%beta
   #Risk<-expit(RR)
   Risk<-pnorm(RR)
   } else Risk<-NULL

   return(Risk)
}


covDF.Nonpar<-function(S, W, Wu){

#########################################
### function to estimate covariate adjusted control distribution
## assuming this is a location scale shift distribution, where mean of marker is linear in covariate

##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
Su<-sort(unique(S))
#S.matrix<-matrix(rep(S,length(Su)),byrow=T,nrow=length(Su))

Dens<-array(NA,dim=c(length(Wu),length(Su)))
for (i in 1:length(Wu)){
 RR<-N.L.E(Su,S[W==Wu[i]])
 RR<-c(0,RR)
 RR<-diff(RR)
 #Dens[i,]<-apply(t(S.matrix==Su)*(W==Wu[i]),2,sum)
 Dens[i,]<-RR/sum(RR)
}
Dens<-as.data.frame(Dens)
colnames(Dens)<-Su
rownames(Dens)<-Wu

#Dens<-Dens/apply(Dens,1,sum)
return(Dens)
}



covDF.Nonpar2<-function(S, W, Wu){

#########################################
### function to estimate covariate adjusted control distribution nonparametrically
## conditional on W but not on Z

##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
SS<-ll<-list()
for (i in 1:length(Wu)){
 SS[[i]]<-sort(S[W==Wu[i]])
 ll[[i]]<-sum(W==Wu[i])
}
return(list(Sval=SS,lval=ll))
}





covDFcc.Nonpar<-function(S, W, Wu, weights){

#########################################
### function to estimate covariate adjusted control distribution
## assuming this is a location scale shift distribution, where mean of marker is linear in covariate

##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
Su<-sort(unique(S))
#S.matrix<-matrix(rep(S,length(Su)),byrow=T,nrow=length(Su))

Dens<-array(NA,dim=c(length(Wu),length(Su)))
for (i in 1:length(Wu)){
 RR<-N.L.E(Su,S[W==Wu[i]])
 RR<-c(0,RR)
 RR<-diff(RR)



 Sw<-S[W==Wu[i]]
 oo<-order(Sw)
 SSo<-Sw[oo]
 WWeights<-tapply((weights[W==Wu[i]])[oo],SSo,sum)
 oo<-match(unique(SSo),Su)
 Sweights<-rep(0,length(Su))
 Sweights[oo]<-WWeights




 RR<-RR*Sweights

 #Dens[i,]<-apply(t(S.matrix==Su)*(W==Wu[i]),2,sum)
 #Dens[i,]<-RR/sum(RR)
 Dens[i,]<-RR
 #######



}
Dens<-as.data.frame(Dens)
colnames(Dens)<-Su
rownames(Dens)<-Wu

#Dens<-Dens/apply(Dens,1,sum)
return(Dens)
}



covDFcc.Nonpar2<-function(S, W, Wu, weights){

#########################################
### function to estimate covariate adjusted control distribution
## nonparametrically, conditional on W only
##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
SS<-ll<-Sweights<-list()
for (i in 1:length(Wu)){
 oo<-order(S[W==Wu[i]])
 SS[[i]]<-(S[W==Wu[i]])[oo]
 ll[[i]]<-sum(W==Wu[i])
 Sweights[[i]]<-(weights[W==Wu[i]])[oo]
}
return(list(Sval=SS,lval=ll,Sweights=Sweights))
}



EM.cc.CPV.Probit.Small.Nonpar.NC<-function(Z,Dens,S1,W,Y,delta,IC,beta,smooth=FALSE,lambda=1,arg=1){

#### smooth: indiator of whether error for location scale family is function of an extra covaraite
#### lambda: lambda for kernel function used for smoothing
#### Ucov: covariate for loc parameter
#### Wcov: covariate for log(scale) parameter

if (arg==1){

ind<-rep(1:sum(delta==0 & IC==0),dim(Dens)[2])  ### for all these subjects in placebo arm and not in IC

YY<-rep(Y[delta==0 & IC==0],dim(Dens)[2])

#S1.res<-rep(S1.out$res,rep(sum(delta==0 & IC==0),sum(delta==1)))


#resWeights<-rep(S1.out$weights,rep(sum(delta==0 & IC==0),sum(delta==1))) ### weights since residuals are not from iid
WW<-rep(W[delta==0 & IC==0],dim(Dens)[2])
ZZ<-rep(Z[delta==0 & IC==0],dim(Dens)[2])



#Loc1<-Ucov[delta==0 & IC==0,]%*%S1.out$eta
#Scale1<-as.matrix(Wcov[delta==0 & IC==0,])%*%S1.out$delta
#
#LLoc1<-rep(Loc1,sum(delta==1))
#SScale1<-rep(Scale1,sum(delta==1))
#SS1<-S1.res*exp(SScale1)+LLoc1

SS1<-as.numeric(rep(colnames(Dens),rep(sum(delta==0 & IC==0), dim(Dens)[2])))





 oo<-match(W[delta==0 & IC==0],rownames(Dens))
 resWeights<-as.numeric(as.matrix(Dens[oo,]))

   # if (smooth==TRUE) {
#        WW.val<-rep(W[delta==1],rep(sum(delta==0 & IC==0),sum(delta==1)))
#        smoothWeights<-exp(-0.5*((WW-WW.val)/lambda)^2)
#        denom.smoothWeights<-tapply(smoothWeights,ind,sum)
#        denom.smoothWeights<-rep(denom.smoothWeights,sum(delta==1))
#        smoothWeights<-smoothWeights/denom.smoothWeights
#
#    }
#
}


beta.old<-beta

#################### arg=2 is not updated yet compared to EM.cc ############################
if (arg==2){

n.C<-sum(Z==1 | IC==1) ### number of complete data
n.R<-sum(IC==0 & Z==0) ### number of reduced data

ind<-rep(1:sum(IC==0 & Z==0),sum(Z==1 | IC==1))  ### for all these subjects in placebo arm and not in IC

ind.Z<-rep(1:n.C,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1))) ### indicator for residual

YY<-rep(Y[IC==0 & Z==0],sum(Z==1 | IC==1))

oo<-IC==1
S1.res.IC<-(S1[oo]-S1.out$eta[1]-S1.out$eta[2]*W[oo])/exp(S1.out$delta)


S1.res<-rep(c(S1.out$res,S1.res.IC),rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))

WW<-rep(W[IC==0 & Z==0],sum(Z==1 | IC==1))
SS1<-S1.res*(exp(S1.out$delta))+S1.out$eta[1]+S1.out$eta[2]*WW


g.old<-rep(1/n.C,n.R)

}
iter<-0
repeat{

### probability of delta=1 conditional on S1,W

### if random sample S from Z=1
#P.Z1<-mean(delta[Z==1 & Y==1])*mean(Y[Z==1])+mean(delta[Z==1 & Y==0])*mean(1-Y[Z==1])

### if random sample S from cases and controls in vaccine arm
P.Z1<-mean(delta[Z==1 & Y==1])*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+mean(delta[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))


#P.Z1<-mean(delta[Z==1])

P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
ff<-ff*resWeights
ff<-ff/PP
if (smooth==TRUE){
    ff<-ff*smoothWeights
}


if (arg==1) {
 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)

 denom<-rep(denom,sum(delta==1))
 weights<-ff/denom
 } else weights<-ff

}  else {

Ljj.old<-tapply(ff*g.old,ind,sum)
Ljj.old<-rep(Ljj.old,sum(Z==1 | IC==1))
g.new<-tapply(ff*g.old/Ljj.old,ind.Z,sum)
g.new<-rep(g.new,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))
g.new<-(1+g.new)/(n.C+n.R)

if (length(ff)>0) {

denom<-tapply(ff*g.new,ind,sum)
denom<-rep(denom,sum(Z==1 | IC==1))
weights<-ff*g.new/denom
} else weights<-ff

}
#YY<-WW<-SS1<-SS2<-weights<-NULL

dataout<-data.frame(Y=c(Y[delta==1 | IC==1],YY),Z=c(Z[delta==1 | IC==1],ZZ),W=c(W[delta==1 | IC==1],WW),
S1=c(S1[delta==1 | IC==1],SS1),weights=c(rep(1,sum(delta==1 | IC==1)),weights))

#date()
#fit<-glm(Y~Z+S1+I(Z*S1)+W+I(Z*W),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
#fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>70 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>70 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}
############33



EM.cc.CPV.Probit.Small.Nonpar2.NC<-function(Z,Sout,S1,W,Y,Wu,delta,IC,beta,smooth=FALSE,lambda=1,arg=1,maxit=200){

#### smooth: indiator of whether error for location scale family is function of an extra covaraite
#### lambda: lambda for kernel function used for smoothing
#### Ucov: covariate for loc parameter
#### Wcov: covariate for log(scale) parameter
Sval<-Sout$Sval
lval<-Sout$lval
if (arg==1){


oo<-match(W[delta==0 & IC==0],Wu)
#ind<-rep(1:sum(delta==0 & IC==0),dim(Dens)[2])  ### for all these subjects in placebo arm and not in IC
ind<-rep(1:sum(delta==0 & IC==0), unlist(lval[oo]))

#YY<-rep(Y[delta==0 & IC==0],dim(Dens)[2])
YY<-rep(Y[delta==0 & IC==0],unlist(lval[oo]))

#S1.res<-rep(S1.out$res,rep(sum(delta==0 & IC==0),sum(delta==1)))


#resWeights<-rep(S1.out$weights,rep(sum(delta==0 & IC==0),sum(delta==1))) ### weights since residuals are not from iid
#WW<-rep(W[delta==0 & IC==0],dim(Dens)[2])
WW<-rep(W[delta==0 & IC==0],unlist(lval[oo]))


#ZZ<-rep(Z[delta==0 & IC==0],dim(Dens)[2])
ZZ<-rep(Z[delta==0 & IC==0],unlist(lval[oo]))


#Loc1<-Ucov[delta==0 & IC==0,]%*%S1.out$eta
#Scale1<-as.matrix(Wcov[delta==0 & IC==0,])%*%S1.out$delta
#
#LLoc1<-rep(Loc1,sum(delta==1))
#SScale1<-rep(Scale1,sum(delta==1))
#SS1<-S1.res*exp(SScale1)+LLoc1

#SS1<-as.numeric(rep(colnames(Dens),rep(sum(delta==0 & IC==0), dim(Dens)[2])))


SS1<-unlist(Sval[oo])


 #oo<-match(W[delta==0 & IC==0],rownames(Dens))
 #resWeight<-Dens[oo,]
 #resWeight<-as.numeric(t(Dens[oo,]))

   # if (smooth==TRUE) {
#        WW.val<-rep(W[delta==1],rep(sum(delta==0 & IC==0),sum(delta==1)))
#        smoothWeights<-exp(-0.5*((WW-WW.val)/lambda)^2)
#        denom.smoothWeights<-tapply(smoothWeights,ind,sum)
#        denom.smoothWeights<-rep(denom.smoothWeights,sum(delta==1))
#        smoothWeights<-smoothWeights/denom.smoothWeights
#
#    }
#
}


beta.old<-beta

#################### arg=2 is not updated yet compared to EM.cc ############################
if (arg==2){

n.C<-sum(Z==1 | IC==1) ### number of complete data
n.R<-sum(IC==0 & Z==0) ### number of reduced data

ind<-rep(1:sum(IC==0 & Z==0),sum(Z==1 | IC==1))  ### for all these subjects in placebo arm and not in IC

ind.Z<-rep(1:n.C,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1))) ### indicator for residual

YY<-rep(Y[IC==0 & Z==0],sum(Z==1 | IC==1))

oo<-IC==1
S1.res.IC<-(S1[oo]-S1.out$eta[1]-S1.out$eta[2]*W[oo])/exp(S1.out$delta)


S1.res<-rep(c(S1.out$res,S1.res.IC),rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))

WW<-rep(W[IC==0 & Z==0],sum(Z==1 | IC==1))
SS1<-S1.res*(exp(S1.out$delta))+S1.out$eta[1]+S1.out$eta[2]*WW


g.old<-rep(1/n.C,n.R)

}
iter<-0
repeat{

### probability of delta=1 conditional on S1,W

### if random sample S from Z=1
#P.Z1<-mean(delta[Z==1 & Y==1])*mean(Y[Z==1])+mean(delta[Z==1 & Y==0])*mean(1-Y[Z==1])

### if random sample S from cases and controls in vaccine arm
P.Z1<-mean(delta[Z==1 & Y==1])*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+mean(delta[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))


#P.Z1<-mean(delta[Z==1])

P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
#ff<-ff*resWeights
ff<-ff/PP
if (smooth==TRUE){
    ff<-ff*smoothWeights
}


if (arg==1) {
 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)

 #denom<-rep(denom,sum(delta==1))
 denom<-denom[ind]
 weights<-as.numeric(ff)/denom
 } else weights<-ff

}  else {

Ljj.old<-tapply(ff*g.old,ind,sum)
Ljj.old<-rep(Ljj.old,sum(Z==1 | IC==1))
g.new<-tapply(ff*g.old/Ljj.old,ind.Z,sum)
g.new<-rep(g.new,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))
g.new<-(1+g.new)/(n.C+n.R)

if (length(ff)>0) {

denom<-tapply(ff*g.new,ind,sum)
denom<-rep(denom,sum(Z==1 | IC==1))
weights<-ff*g.new/denom
} else weights<-ff

}
#YY<-WW<-SS1<-SS2<-weights<-NULL

dataout<-data.frame(Y=c(Y[delta==1 | IC==1],YY),Z=c(Z[delta==1 | IC==1],ZZ),W=c(W[delta==1 | IC==1],WW),
S1=c(S1[delta==1 | IC==1],SS1),weights=c(rep(1,sum(delta==1 | IC==1)),weights))

#date()
#fit<-glm(Y~Z+S1+I(Z*S1)+W+I(Z*W),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
#fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>maxit | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>maxit | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}


########################3 What if we don't utilize the condition that X|W=X|Z,W ##################3




covDFcc.Nonpar2.ZW<-function(S, W, Wu, Z,Zu,weights){

#########################################
### function to estimate covariate adjusted control distribution
## nonparametrically, conditional Z and W
##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
SS<-ll<-Sweights<-list()
ZWu<-expand.grid(Zu,Wu)
for (i in 1:dim(ZWu)[1]){
 oo<-order(S[Z==ZWu[i,1] & W==ZWu[i,2]])
 SS[[i]]<-(S[Z==ZWu[i,1] & W==ZWu[i,2]])[oo]
 ll[[i]]<-sum(Z==ZWu[i,1] & W==ZWu[i,2])
 Sweights[[i]]<-(weights[Z==ZWu[i,1] & W==ZWu[i,2]])[oo]
}
return(list(Sval=SS,lval=ll,Sweights=Sweights))
}




EM.cc.CPV.Probit.Small.Nonpar2.NC.ZW<-function(Z,Sout,S1,W,Y,Zu,Wu,delta,IC,beta,smooth=FALSE,lambda=1,arg=1){

#### smooth: indiator of whether error for location scale family is function of an extra covaraite
#### lambda: lambda for kernel function used for smoothing
#### Ucov: covariate for loc parameter
#### Wcov: covariate for log(scale) parameter
Sval<-Sout$Sval
lval<-Sout$lval

ZWu<-expand.grid(Zu,Wu)
ZWu<-ZWu[,1]*10+ZWu[,2]

ZW<-Z*10+W
if (arg==1){


oo<-match(ZW[delta==0 & IC==0],ZWu)
#ind<-rep(1:sum(delta==0 & IC==0),dim(Dens)[2])  ### for all these subjects in placebo arm and not in IC
ind<-rep(1:sum(delta==0 & IC==0), unlist(lval[oo]))

#YY<-rep(Y[delta==0 & IC==0],dim(Dens)[2])
YY<-rep(Y[delta==0 & IC==0],unlist(lval[oo]))

#S1.res<-rep(S1.out$res,rep(sum(delta==0 & IC==0),sum(delta==1)))


#resWeights<-rep(S1.out$weights,rep(sum(delta==0 & IC==0),sum(delta==1))) ### weights since residuals are not from iid
#WW<-rep(W[delta==0 & IC==0],dim(Dens)[2])
WW<-rep(W[delta==0 & IC==0],unlist(lval[oo]))


#ZZ<-rep(Z[delta==0 & IC==0],dim(Dens)[2])
ZZ<-rep(Z[delta==0 & IC==0],unlist(lval[oo]))


#Loc1<-Ucov[delta==0 & IC==0,]%*%S1.out$eta
#Scale1<-as.matrix(Wcov[delta==0 & IC==0,])%*%S1.out$delta
#
#LLoc1<-rep(Loc1,sum(delta==1))
#SScale1<-rep(Scale1,sum(delta==1))
#SS1<-S1.res*exp(SScale1)+LLoc1

#SS1<-as.numeric(rep(colnames(Dens),rep(sum(delta==0 & IC==0), dim(Dens)[2])))


SS1<-unlist(Sval[oo])


 #oo<-match(W[delta==0 & IC==0],rownames(Dens))
 #resWeight<-Dens[oo,]
 #resWeight<-as.numeric(t(Dens[oo,]))

   # if (smooth==TRUE) {
#        WW.val<-rep(W[delta==1],rep(sum(delta==0 & IC==0),sum(delta==1)))
#        smoothWeights<-exp(-0.5*((WW-WW.val)/lambda)^2)
#        denom.smoothWeights<-tapply(smoothWeights,ind,sum)
#        denom.smoothWeights<-rep(denom.smoothWeights,sum(delta==1))
#        smoothWeights<-smoothWeights/denom.smoothWeights
#
#    }
#
}


beta.old<-beta

#################### arg=2 is not updated yet compared to EM.cc ############################
if (arg==2){

n.C<-sum(Z==1 | IC==1) ### number of complete data
n.R<-sum(IC==0 & Z==0) ### number of reduced data

ind<-rep(1:sum(IC==0 & Z==0),sum(Z==1 | IC==1))  ### for all these subjects in placebo arm and not in IC

ind.Z<-rep(1:n.C,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1))) ### indicator for residual

YY<-rep(Y[IC==0 & Z==0],sum(Z==1 | IC==1))

oo<-IC==1
S1.res.IC<-(S1[oo]-S1.out$eta[1]-S1.out$eta[2]*W[oo])/exp(S1.out$delta)


S1.res<-rep(c(S1.out$res,S1.res.IC),rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))

WW<-rep(W[IC==0 & Z==0],sum(Z==1 | IC==1))
SS1<-S1.res*(exp(S1.out$delta))+S1.out$eta[1]+S1.out$eta[2]*WW


g.old<-rep(1/n.C,n.R)

}
iter<-0
repeat{

### probability of delta=1 conditional on S1,W

### if random sample S from Z=1
#P.Z1<-mean(delta[Z==1 & Y==1])*mean(Y[Z==1])+mean(delta[Z==1 & Y==0])*mean(1-Y[Z==1])

### if random sample S from cases and controls in vaccine arm
#P.Z1<-mean(delta[Z==1 & Y==1])*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+mean(delta[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))
#
#
##P.Z1<-mean(delta[Z==1])
#
#P.Z0<-mean(delta[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
#PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

###############

P.Y1Z1<-mean(delta[Z==1 & Y==1])
P.Y1Z0<-mean(delta[Z==0 & Y==1])
P.Y0Z0<-mean(delta[Z==0 & Y==0])
P.Y0Z1<-mean(delta[Z==1 & Y==0])

PP<-(ZZ==1)*(P.Y1Z1*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+P.Y0Z1*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)))+
(ZZ==0)*(P.Y1Z0*R.probit(Z=0,S1=SS1,W=WW,beta=beta.old)+P.Y0Z0*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old)))



ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
#ff<-ff*resWeights
ff<-ff/PP
if (smooth==TRUE){
    ff<-ff*smoothWeights
}


if (arg==1) {
 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)

 #denom<-rep(denom,sum(delta==1))
 denom<-denom[ind]
 weights<-as.numeric(ff)/denom
 } else weights<-ff

}  else {

Ljj.old<-tapply(ff*g.old,ind,sum)
Ljj.old<-rep(Ljj.old,sum(Z==1 | IC==1))
g.new<-tapply(ff*g.old/Ljj.old,ind.Z,sum)
g.new<-rep(g.new,rep(sum(IC==0 & Z==0),sum(Z==1 | IC==1)))
g.new<-(1+g.new)/(n.C+n.R)

if (length(ff)>0) {

denom<-tapply(ff*g.new,ind,sum)
denom<-rep(denom,sum(Z==1 | IC==1))
weights<-ff*g.new/denom
} else weights<-ff

}
#YY<-WW<-SS1<-SS2<-weights<-NULL

dataout<-data.frame(Y=c(Y[delta==1 | IC==1],YY),Z=c(Z[delta==1 | IC==1],ZZ),W=c(W[delta==1 | IC==1],WW),
S1=c(S1[delta==1 | IC==1],SS1),weights=c(rep(1,sum(delta==1 | IC==1)),weights))

#date()
#fit<-glm(Y~Z+S1+I(Z*S1)+W+I(Z*W),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
#fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)
fit<-glm(Y~Z+S1+I(Z*S1),weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>200 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>200 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}



##### allows case-control sampling of S among vaccinees also

Ssub.cc.Alter<-function(datain,ratio.v,ratio.p){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


if (ratio.p==-1){
    p.IC.1=1
} else {
p.IC.1<-ratio.p
}

p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)

IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==0 | p.IC==1) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


if (ratio.v==-1){
 p.CC.1=1
} else {
p.CC.1<-ratio.v
}
p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

#### do an extra sampling just to keep the seed same #######
if (p.CC==0 | p.CC==1) rbinom(N,1,0.1)

weights<-ifelse(Z==0,0,1)
weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
IC.vac<-IC.vac | (Z==1 & Y==1)


################# select subcohort for IC ######

datain$IC<-IC
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}




##### allows case-control sampling of S among vaccinees also

Ssub.cc.Alter.Full<-function(datain,ratio.v,ratio.p,ratio.vc=1){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


if (ratio.p==-1){
    p.IC.1=1
} else {
p.IC.1<-ratio.p
}

p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)

IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==1 | p.IC==0) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


if (ratio.v==-1){
 p.CC.1=1
} else {
p.CC.1<-ratio.v
}
p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)


#### do an extra sampling just to keep the seed same #######
if (p.CC==1 | p.CC==0) rbinom(N,1,0.1)


#############

if (ratio.vc==1)
{
rbinom(N,1,0.1)
IC.vac.case<-(Z==1) & (Y==1)
}
else {
  IC.vac.case<-((Z==1) & (Y==1))*rbinom(N,1,ratio.vc)
}


weights<-ifelse(Z==0,0,1)
weights<-ifelse(Z==1 & Y==1, 1/(sum(IC.vac.case)/sum(Z==1 & Y==1)),weights)
weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
IC.vac<-IC.vac | IC.vac.case


################# select subcohort for IC ######

datain$IC<-IC
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}








##### allows case-control sampling of S among vaccinees also

### allow for sampling of cases among placebo
Ssub.cc.Alter.Full.Obs<-function(datain,ratio.v,ratio.p,ratio.vc=1,ratio.pc=1){

##### sample subcohort #################
N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


if (ratio.p==-1){
    p.IC.1=1
} else {
p.IC.1<-ratio.p
}

p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)

IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
if (p.IC==1 | p.IC==0) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


if (ratio.v==-1){
 p.CC.1=1
} else {
p.CC.1<-ratio.v
}
p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)


#### do an extra sampling just to keep the seed same #######
if (p.CC==1 | p.CC==0) rbinom(N,1,0.1)


#############

if (ratio.vc==1)
{
rbinom(N,1,0.1)
IC.vac.case<-(Z==1) & (Y==1)
}
else {
  IC.vac.case<-((Z==1) & (Y==1))*rbinom(N,1,ratio.vc)
}

if (ratio.pc==1)
{
rbinom(N,1,0.1)
IC.plc.case<-(Z==0) & (Y==1)
}
else {
  IC.plc.case<-((Z==0) & (Y==1))*rbinom(N,1,ratio.pc)
}

#weights<-ifelse(Z==0,0,1)

weights<-ifelse(Z==0 & Y==1, 1/(sum(IC.plc.case)/sum(Z==0 & Y==1)),0)
weights<-ifelse(Z==0 & Y==0, 1/(sum(IC)/sum(Z==0 & Y==0)),weights)
weights<-ifelse(Z==1 & Y==1, 1/(sum(IC.vac.case)/sum(Z==1 & Y==1)),weights)
weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
IC.vac<-IC.vac | IC.vac.case


################# select subcohort for IC ######

datain$IC<-IC | IC.plc.case
datain$IC.vac<-IC.vac
datain$Y<-Y
datain$S1<-ifelse(datain$IC | datain$IC.vac, datain$S1,NA)
datain$S2<-ifelse(datain$IC | datain$IC.vac, datain$S2,NA)
datain$weights<-weights
return(dataout=datain)
}




















rmulti<-function(N,ratio.w,ratio){
 if (ratio.w<ratio) return(matrix(NA,3,N))
 out1=rbinom(N,1,ratio.w)
 if (ratio.w==0) {
   rbinom(N,1,.1)
   rbinom(N,1,.1)
   out2=rep(0,N)
 } else if (ratio.w==1){
    rbinom(N,1,.1)
    out2=rbinom(N,1,ratio/ratio.w)*out1
    if (ratio==ratio.w | ratio==0) rbinom(N,1,.1)
  } else {
    out2=rbinom(N,1,ratio/ratio.w)*out1
    if (ratio==ratio.w | ratio==0) rbinom(N,1,.1)
  }
  out=rbind(out1-out2,out2,1-out1)
  return(out)
 }


##### allows case-control sampling of S among vaccinees also, and allow subsampling of W, but S has to be nested within W

Ssub.cc.Alter.Full.SubW.old<-function(datain,ratio.v.w, ratio.p.w, ratio.vc.w, ratio.pc.w, ratio.v,ratio.p,ratio.vc=1){

### ratio.vc.w is probability of sampling W among vaccinated cases
### ratio.vc is probability of sampling S among vaccinated cases
### ratio.v.w is probability of sampling W among vaccinated controls
### ratio.v is probability of sampling S among vaccinated controls
### ratio.pc.w is probability of sampling W among placebo cases
### ratio.p.w is probability of sampling W among placebo controls
### ratio.p is probability of sampling S among placebo controls

##### sample subcohort #################

#ratio.vc.w=1
#ratio.vc=1
#ratio.v.w=0.3
#ratio.v=0.3
#ratio.pc.w=1
ratio.pc=0
#ratio.p.w=0.3
#ratio.p=0.3

N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


#if (ratio.p==-1){
#    p.IC.1=1
#} else {
#p.IC.1<-ratio.p
#}
#
#p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)
#

IC.Y0Z0<-t(((Z==0)*(Y==0))*t(rmultinom(N,1,c(ratio.p.w-ratio.p,ratio.p,1-ratio.p.w))))




#IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
#if (p.IC==1 | p.IC==0) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


#if (ratio.v==-1){
# p.CC.1=1
#} else {
#p.CC.1<-ratio.v
#}
#p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
#
#IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

IC.Y1Z0<-t(((Z==0)*(Y==1))*t(rmultinom(N,1,c(ratio.pc.w-ratio.pc,ratio.pc,1-ratio.pc.w))))

IC.Y0Z1<-t(((Z==1)*(Y==0))*t(rmultinom(N,1,c(ratio.v.w-ratio.v,ratio.v,1-ratio.v.w))))

IC.Y1Z1<-t(((Z==1)*(Y==1))*t(rmultinom(N,1,c(ratio.vc.w-ratio.vc,ratio.vc,1-ratio.vc.w))))

#### do an extra sampling just to keep the seed same #######
#if (p.CC==1 | p.CC==0) rbinom(N,1,0.1)


#############

#if (ratio.vc==1)
#{
#rbinom(N,1,0.1)
#IC.vac.case<-(Z==1) & (Y==1)
#}
#else {
#  IC.vac.case<-((Z==1) & (Y==1))*rbinom(N,1,ratio.vc)
#}

delta1=(IC.Y0Z0[1,]+IC.Y0Z0[2,]>0)|(IC.Y0Z1[1,]+IC.Y0Z1[2,]>0)|(IC.Y1Z0[1,]+IC.Y1Z0[2,]>0)|(IC.Y1Z1[1,]+IC.Y1Z1[2,]>0)
delta2=(IC.Y0Z0[2,]>0)|(IC.Y0Z1[2,]>0)|(IC.Y1Z0[2,]>0)|(IC.Y1Z1[2,])


#weights<-ifelse(Z==0,0,1)
#weights<-ifelse(Z==1 & Y==1, 1/(sum(IC.vac.case)/sum(Z==1 & Y==1)),weights)
#weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
#IC.vac<-IC.vac | IC.vac.case


################# select subcohort for IC ######

datain$delta1=delta1
datain$delta2=delta2
datain$Y<-Y
datain$W<-ifelse(datain$delta1, datain$W, NA)
datain$S1<-ifelse(datain$delta2, datain$S1,NA)
datain$S2<-ifelse(datain$delta2, datain$S2,NA)
#datain$weights<-weights
return(dataout=datain)
}



##### allows case-control sampling of S among vaccinees also, and allow subsampling of W, but S has to be nested within W
### this code do extra sampling to keep the right seed
Ssub.cc.Alter.Full.SubW<-function(datain,ratio.v.w, ratio.p.w, ratio.vc.w, ratio.pc.w, ratio.v,ratio.p,ratio.vc=1){

### ratio.vc.w is probability of sampling W among vaccinated cases
### ratio.vc is probability of sampling S among vaccinated cases
### ratio.v.w is probability of sampling W among vaccinated controls
### ratio.v is probability of sampling S among vaccinated controls
### ratio.pc.w is probability of sampling W among placebo cases
### ratio.p.w is probability of sampling W among placebo controls
### ratio.p is probability of sampling S among placebo controls

##### sample subcohort #################

#ratio.vc.w=1
#ratio.vc=1
#ratio.v.w=0.3
#ratio.v=0.3
#ratio.pc.w=1
ratio.pc=0
#ratio.p.w=0.3
#ratio.p=0.3

N<-dim(datain)[1]

R<-datain$R
Z<-datain$Z
#S1<-datain$S1
#S2<-datain$S2
#W<-datain$W

Y<-rbinom(N,1,R)


#if (ratio.p==-1){
#    p.IC.1=1
#} else {
#p.IC.1<-ratio.p
#}
#
#p.IC<-ifelse(p.IC.1>1, 1, p.IC.1)
#

IC.Y0Z0<-t(((Z==0)*(Y==0))*t(rmulti(N,ratio.p.w,ratio.p)))



#IC<-((Z==0) & (Y==0))*rbinom(N,1,p.IC)

#### do an extra sampling just to keep the seed same  ####
#if (p.IC==1 | p.IC==0) rbinom(N,1,0.1)
### assume the mean # of controls sampled is 5 times that of cases ####
### could be modified if want the exact # of controls = 5* # of  cases


#if (ratio.v==-1){
# p.CC.1=1
#} else {
#p.CC.1<-ratio.v
#}
#p.CC<-ifelse(p.CC.1>1, 1, p.CC.1)
#
#IC.vac<-((Z==1) & (Y==0))*rbinom(N,1,p.CC)

IC.Y1Z0<-t(((Z==0)*(Y==1))*t(rmulti(N,ratio.pc.w,ratio.pc)))

IC.Y0Z1<-t(((Z==1)*(Y==0))*t(rmulti(N,ratio.v.w,ratio.v)))

IC.Y1Z1<-t(((Z==1)*(Y==1))*t(rmulti(N,ratio.vc.w,ratio.vc)))

#### do an extra sampling just to keep the seed same #######
#if (p.CC==1 | p.CC==0) rbinom(N,1,0.1)


#############

#if (ratio.vc==1)
#{
#rbinom(N,1,0.1)
#IC.vac.case<-(Z==1) & (Y==1)
#}
#else {
#  IC.vac.case<-((Z==1) & (Y==1))*rbinom(N,1,ratio.vc)
#}

delta1=(IC.Y0Z0[1,]+IC.Y0Z0[2,]>0)|(IC.Y0Z1[1,]+IC.Y0Z1[2,]>0)|(IC.Y1Z0[1,]+IC.Y1Z0[2,]>0)|(IC.Y1Z1[1,]+IC.Y1Z1[2,]>0)
delta2=(IC.Y0Z0[2,]>0)|(IC.Y0Z1[2,]>0)|(IC.Y1Z0[2,]>0)|(IC.Y1Z1[2,])


#weights<-ifelse(Z==0,0,1)
#weights<-ifelse(Z==1 & Y==1, 1/(sum(IC.vac.case)/sum(Z==1 & Y==1)),weights)
#weights<-ifelse(Z==1 & Y==0, 1/(sum(IC.vac)/sum(Z==1 & Y==0)),weights)
#IC.vac<-IC.vac | IC.vac.case


################# select subcohort for IC ######

datain$delta1=delta1
datain$delta2=delta2
datain$Y<-Y
datain$W<-ifelse(datain$delta1, datain$W, NA)
datain$S1<-ifelse(datain$delta2, datain$S1,NA)
datain$S2<-ifelse(datain$delta2, datain$S2,NA)
#datain$weights<-weights
return(dataout=datain)
}

#Z=Z;Sout=Sout.NC.All;S1=S1;W=W;Y=Y;Wu=Wu;delta1=delta1;delta2=delta2;p.delta2=p.delta2;pd2.YZW=pd2.YZW;beta=gamma

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.old<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,p.delta2,pd2.YZW,beta){

    Sval<-Sout$Sval
    lval<-Sout$lval

    oo.2<-match(W[delta1==1 & delta2==0],Wu)

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(WW.2,Wu)
    #############
    oo.1<-rep(Wu,unlist(lval))

    ind.1<-rep(1:sum(delta1==0), sum(unlist(lval)))
    YY.1<-rep(Y[delta1==0],sum(unlist(lval)))
    WW.1<-rep(Wu,sum(delta1==0)*unlist(lval))
    ZZ.1<-rep(Z[delta1==0],sum(unlist(lval)))
    SS1.1<-rep(unlist(Sval),sum(delta1==0))


    X<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))

    oo.11<-matrix(WW.1,Wu)


    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W

    ### if random sample S from cases and controls in vaccine arm

    ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####

    ### here is a possible bug #### check, oo.2 is of length smaller than SS1.2
    P.Z1<-pd2.YZW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old)+pd2.YZW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old))

    P.Z0<-pd2.YZW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W

     ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-pd2.YZW[2,2,oo.11]*R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old)+pd2.YZW[1,2,oo.11]*(1-R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old))

    P.Z0<-pd2.YZW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

    ### probability of delta2==1 conditional on delta1==1 and W

    PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff/PP*PP.delta2



     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(X,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,p.delta2,pd2.YZW,beta){

#Sout=Sout.NC.All;Wu=Wu=c(1,2,3,4);beta=gamma



   Wu.weights<-rep(NA,length(Wu))
   for (i in 1:length(Wu)){
        Wu.weights[i]=mean(W[delta1==1]==Wu[i])
   }

    Sval<-Sout$Sval
    lval<-Sout$lval

    oo.2<-match(W[delta1==1 & delta2==0],Wu)

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(WW.2,Wu)
    #############
    oo.1<-rep(Wu,unlist(lval))


    ind.1<-rep(1:sum(delta1==0), sum(unlist(lval)))
    YY.1<-rep(Y[delta1==0],sum(unlist(lval)))
    WW.1<-rep(Wu,sum(delta1==0)*unlist(lval))
    ZZ.1<-rep(Z[delta1==0],sum(unlist(lval)))
    #### note the following about SS1.1 has corrected a bug pointed out by Krisz on 4/19
    SS1.1<-rep(unlist(Sval),each=sum(delta1==0))
    Sweights.1<-rep(Wu.weights,sum(delta1==0)*unlist(lval))




    X<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))

    oo.11<-matrix(WW.1,Wu)


    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W

    ### if random sample S from cases and controls in vaccine arm

    ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####

    ### here is a possible bug #### check, oo.2 is of length smaller than SS1.2
    P.Z1<-pd2.YZW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old)+pd2.YZW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old))

    P.Z0<-pd2.YZW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W

     ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-pd2.YZW[2,2,oo.11]*R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old)+pd2.YZW[1,2,oo.11]*(1-R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old))

    P.Z0<-pd2.YZW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

    ### probability of delta2==1 conditional on delta1==1 and W

    PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff/PP*PP.delta2*Sweights.1



     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(X,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}

### this function is a special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z)  ######################
EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.special.old<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,p.delta2,beta){

    Sval<-Sout$Sval
    lval<-Sout$lval

    oo.2<-match(W[delta1==1 & delta2==0],Wu)

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    #############
    oo.1<-rep(Wu,unlist(lval))

    ind.1<-rep(1:sum(delta1==0), sum(unlist(lval)))
    YY.1<-rep(Y[delta1==0],sum(unlist(lval)))
    WW.1<-rep(Wu,sum(delta1==0)*unlist(lval))
    ZZ.1<-rep(Z[delta1==0],sum(unlist(lval)))
    SS1.1<-rep(unlist(Sval),sum(delta1==0))


    X<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))



    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W

    ### if random sample S from cases and controls in vaccine arm

    ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old))

    P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1.2,W=WW.2,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W

     ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old))

    P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1.1,W=WW.1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

    ### probability of delta2==1 conditional on delta1==1 and W

    PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff/PP*PP.delta2



     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(X,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.special<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,p.delta2,beta){

   Wu.weights<-rep(NA,length(Wu))
   for (i in 1:length(Wu)){
        Wu.weights[i]=mean(W[delta1==1]==Wu[i])
   }

    Sval<-Sout$Sval
    lval<-Sout$lval

    oo.2<-match(W[delta1==1 & delta2==0],Wu)

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    #############
    oo.1<-rep(Wu,unlist(lval))

    ind.1<-rep(1:sum(delta1==0), sum(unlist(lval)))
    YY.1<-rep(Y[delta1==0],sum(unlist(lval)))
    WW.1<-rep(Wu,sum(delta1==0)*unlist(lval))
    ZZ.1<-rep(Z[delta1==0],sum(unlist(lval)))
    SS1.1<-rep(unlist(Sval),each=sum(delta1==0))
    Sweights.1<-rep(Wu.weights,sum(delta1==0)*unlist(lval))


    X<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))



    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W

    ### if random sample S from cases and controls in vaccine arm

    ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1.2,W=WW.2,beta=beta.old))

    P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1.2,W=WW.2,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W

     ### for the special case where P(delta1=delta2=1|Y,Z,W)=P(delta1=delta2=1|Y,Z) ####
    P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit(Z=1,S1=SS1.1,W=WW.1,beta=beta.old))

    P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit(Z=0,S1=SS1.1,W=WW.1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

    ### probability of delta2==1 conditional on delta1==1 and W

    PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff/PP*PP.delta2*Sweights.1




     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(X,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.L.special<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,p.delta2,beta){

Sval<-Sout$Sval
lval<-Sout$lval

oo.2<-match(W[delta1==1 & delta2==0],Wu)

ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval[oo.2]))
YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval[oo.2]))
WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval[oo.2]))
ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval[oo.2]))
SS1.2<-unlist(Sval[oo.2])

#############
oo.1<-rep(Wu,unlist(lval))

ind.1<-rep(1:sum(delta1==0), sum(unlist(lval)))
YY.1<-rep(Y[delta1==0],sum(unlist(lval)))
WW.1<-rep(Wu,sum(delta1==0)*unlist(lval))
ZZ.1<-rep(Z[delta1==0],sum(unlist(lval)))
SS1.1<-rep(unlist(Sval),sum(delta1==0))





beta.old<-beta

iter=0
repeat{

### probability of delta2=1 conditional on S1,W

### if random sample S from cases and controls in vaccine arm
P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit.L(Z=1,S1=SS1.2,W=WW.2,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit.L(Z=1,S1=SS1.2,W=WW.2,beta=beta.old))

P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit.L(Z=0,S1=SS1.2,W=WW.2,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


ff<-R.probit.L(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
ff<-ifelse(YY.2==1,ff,1-ff)
ff<-ff/PP



 if (length(ff)>0){
 denom<-tapply(ff,ind.2,sum)
 denom<-denom[ind.2]
 weights.2<-as.numeric(ff)/denom
 } else weights.2<-ff


### probability of delta2=1 conditional on S1,W

### if random sample S from cases and controls in vaccine arm
P.Z1<-mean(delta2[Z==1 & Y==1])*R.probit.L(Z=1,S1=SS1.1,W=WW.1,beta=beta.old)+mean(delta2[Z==1 & Y==0])*(1-R.probit.L(Z=1,S1=SS1.1,W=WW.1,beta=beta.old))

P.Z0<-mean(delta2[Z==0 & Y==0])*(1-R.probit.L(Z=0,S1=SS1.1,W=WW.1,beta=beta.old))
PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

### probability of delta2==1 conditional on delta1==1 and W

PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))

ff<-R.probit.L(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
ff<-ifelse(YY.1==1,ff,1-ff)
ff<-ff/PP*PP.delta2



 if (length(ff)>0){
 denom<-tapply(ff,ind.1,sum)
 denom<-denom[ind.1]
 weights.1<-as.numeric(ff)/denom
 } else weights.1<-ff

dataout<-data.frame(Y=c(Y[delta2==1],YY.2,YY.1),Z=c(Z[delta2==1],ZZ.2,ZZ.1),W=c(W[delta2==1],WW.2,WW.1),
S1=c(S1[delta2==1],SS1.2,SS1.1),weights=c(rep(1,sum(delta2==1)),weights.2,weights.1))

fit<-glm(Y~Z+S1+I(Z*S1)+W,weights=weights,data=dataout,family=binomial(link=probit),start=beta.old)

#date()
#date()

out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

#print(iter);
#print(beta.new);
if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new

#print(iter)
if (iter>200 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>200 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}




  calPd1.ZY.sub<-function(Y,Z,delta1,y,z){

    #### Note  #####
    #### Captial Z,Y,delta1 are data input for calculating the probability of delta1=1 given Y=y, Z=z
     #

     #sum(delta1==1 & Y==y & Z==z)/sum(Y==y & Z==z)
     apply(cbind(z,y),1,function(ci)
       sum(delta1==1 & Y==ci[2] & Z==ci[1])/sum(Y==ci[2] & Z==ci[1])
     )

    }

    ### compute P(delta2=1|delta1=1,Y,Z,W) ###
    calPd2.d1YZW.sub<-function(W,Y,Z,delta1,delta2,w,y,z){

    ### W,Z,Y,delta2,delta1 are data input for calculating the probability of delta2==1 given W=w, Z=z,Y=y, and delta1=1

      #sum(delta2==1 & delta1==1 & Y==y & Z==z & W==w)/sum(delta1==1 & Y==y & Z==z & W==w)

      apply(cbind(w,z,y),1,function(ci)
       sum(delta2==1 & delta1==1 & Y==ci[3] & Z==ci[2] & W==ci[1])/sum(delta1==1 & Y==ci[3] & Z==ci[2] & W==ci[1])
     )

    }

EM.cc.CPV.Probit.Small.Nonpar2.NC.Weighted<-function(Z,Sout,S1,W,Y,Wu,delta1,delta2,beta,smooth=FALSE,lambda=1,arg=1){



        pi.Y0Z0=sum(delta1==1 & Y==0 & Z==0)/sum(Y==0 & Z==0)
        pi.Y1Z0=sum(delta1==1 & Y==1 & Z==0)/sum(Y==1 & Z==0)
        pi.Y0Z1=sum(delta1==1 & Y==0 & Z==1)/sum(Y==0 & Z==1)
        pi.Y1Z1=sum(delta1==1 & Y==1 & Z==1)/sum(Y==1 & Z==1)

        pi.YZ<-c(pi.Y0Z0,pi.Y1Z0,pi.Y0Z1,pi.Y1Z1)




        pi.Y1Z1Wu<-calPd1.ZY.sub(Y,Z,delta1,1,1)*calPd2.d1YZW.sub(W,Y,Z,delta1,delta2,Wu,1,1)
        pi.Y0Z1Wu<-calPd1.ZY.sub(Y,Z,delta1,0,1)*calPd2.d1YZW.sub(W,Y,Z,delta1,delta2,Wu,0,1)
        pi.Y1Z0Wu<-calPd1.ZY.sub(Y,Z,delta1,1,0)*calPd2.d1YZW.sub(W,Y,Z,delta1,delta2,Wu,1,0)
        pi.Y0Z0Wu<-calPd1.ZY.sub(Y,Z,delta1,0,0)*calPd2.d1YZW.sub(W,Y,Z,delta1,delta2,Wu,0,0)

         delta=delta2[delta1]
         Z=Z[delta1]
         Y=Y[delta1]
         W=W[delta1]
         S1=S1[delta1]







    Sval<-Sout$Sval
    lval<-Sout$lval
    if (arg==1){


    oo<-match(W[delta==0],Wu)

    ind<-rep(1:sum(delta==0), unlist(lval[oo]))


    YY<-rep(Y[delta==0],unlist(lval[oo]))


    WW<-rep(W[delta==0],unlist(lval[oo]))


    ZZ<-rep(Z[delta==0],unlist(lval[oo]))




    SS1<-unlist(Sval[oo])



    }


    beta.old<-beta



    oo.w<-match(WW,Wu)
    pi.Y0Z0WW=pi.Y0Z0Wu[oo.w]
    pi.Y1Z0WW=pi.Y1Z0Wu[oo.w]
    pi.Y0Z1WW=pi.Y0Z1Wu[oo.w]
    pi.Y1Z1WW=pi.Y1Z1Wu[oo.w]

    YZu<-expand.grid(c(0,1),c(0,1))
    YZu=YZu[,1]*10+YZu[,2]
    YYZZ<-c(Y[delta==1],YY)*10+c(Z[delta==1],ZZ)

    oo.yz<-match(YYZZ,YZu)

    pi.YYZZ<-pi.YZ[oo.yz]

    X<-as.matrix(cbind(1,c(Z[delta==1],ZZ),c(S1[delta==1],SS1),c(Z[delta==1],ZZ)*c(S1[delta==1],SS1)))
    Y.long=as.vector(c(Y[delta==1],YY))

    iter<-0
    repeat{

    ### probability of delta=1 conditional on S1,W


    ### if random sample S from cases and controls in vaccine arm
    P.Z1<-pi.Y1Z1WW*R.probit(Z=1,S1=SS1,W=WW,beta=beta.old)+pi.Y0Z1WW*(1-R.probit(Z=1,S1=SS1,W=WW,beta=beta.old))


    P.Z0<-pi.Y0Z0WW*(1-R.probit(Z=0,S1=SS1,W=WW,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
    ff<-ifelse(YY==1,ff,1-ff)

    ff<-ff/PP
    if (smooth==TRUE){
        ff<-ff*smoothWeights
    }


    if (arg==1) {
     if (length(ff)>0){
     denom<-tapply(ff,ind,sum)

     denom<-denom[ind]
     weights<-as.numeric(ff)/denom
     } else weights<-ff

    }


    fit=try.error(glm.fit(X,Y.long,weights=c(rep(1,sum(delta==1)),weights)/pi.YYZZ,family=binomial(link=probit),start=beta.old))

    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3


    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}



EM.cc.CPV.Probit.Small.Nonpar2.NC.Weighted.X.Short<-function(Z,Sout,S1,W,Y,X,Wu,Xu,delta1,delta2,beta,special=0,smooth=FALSE,lambda=1,arg=1,kk.RV144=NULL){

 #Sout=Sout.NC.All.X;beta=gamma;special=0;smooth=FALSE;lambda=1;arg=1;kk.RV144=NULL

        #pi.Y0Z0=sum(delta1==1 & Y==0 & Z==0)/sum(Y==0 & Z==0)
        #        pi.Y1Z0=sum(delta1==1 & Y==1 & Z==0)/sum(Y==1 & Z==0)
        #        pi.Y0Z1=sum(delta1==1 & Y==0 & Z==1)/sum(Y==0 & Z==1)
        #        pi.Y1Z1=sum(delta1==1 & Y==1 & Z==1)/sum(Y==1 & Z==1)
        #
        #        pi.YZ<-c(pi.Y0Z0,pi.Y1Z0,pi.Y0Z1,pi.Y1Z1)



        pi.Y0Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,Xu)
        pi.Y1Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,Xu)
        pi.Y0Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,Xu)
        pi.Y1Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,Xu)


        pi.XYZ<-c(pi.Y0Z0Xu,pi.Y1Z0Xu,pi.Y0Z1Xu,pi.Y1Z1Xu)


        WXu.ep<-expand.grid(Wu,Xu)

        if (special==0){

            pi.Y1Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,1)
            pi.Y0Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,1)
            pi.Y1Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,0)
            pi.Y0Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,0)
        } else {
       ### for the special case where p(delta2=1|Y,Z,W,S,X)=p(delta2=1|Y,Z,X)  & p(delta2=1|Y,Z=0,X)=0 ###
           kk.special=table(delta2,Y,Z,X)
           pd2.YZXW<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
           for (j in 1:length(Xu)){
           for (k in 1:length(Wu)){
             pd2.YZXW[,,(j-1)*length(Wu)+k]=kk.RV144[2,,,j]/apply(kk.special[,,,j],c(2,3),sum)
            }
            }

       pd2.YZXW[,1,]<-0

       pi.Y0Z0WXu<-pd2.YZXW[1,1,]
       pi.Y1Z0WXu<-pd2.YZXW[2,1,]
       pi.Y0Z1WXu<-pd2.YZXW[1,2,]
       pi.Y1Z1WXu<-pd2.YZXW[2,2,]
       }

        WXu<-WXu.ep[,1]+WXu.ep[,2]*10

         delta=delta2[delta1]
         Z=Z[delta1]
         Y=Y[delta1]
         W=W[delta1]
         X=X[delta1]
         S1=S1[delta1]


         WX<-W+X*10

    Sval<-Sout$Sval
    lval<-Sout$lval


    if (arg==1){


    oo<-match(WX[delta==0],WXu)

    ind<-rep(1:sum(delta==0), unlist(lval[oo]))


    YY<-rep(Y[delta==0],unlist(lval[oo]))


    WW<-rep(W[delta==0],unlist(lval[oo]))


    ZZ<-rep(Z[delta==0],unlist(lval[oo]))

    XX<-rep(X[delta==0],unlist(lval[oo]))


    SS1<-unlist(Sval[oo])

    WWXX<-WW+XX*10

    }


    beta.old<-beta



    oo.wx<-match(WWXX,WXu)
    pi.Y0Z0WWXX=pi.Y0Z0WXu[oo.wx]
    pi.Y1Z0WWXX=pi.Y1Z0WXu[oo.wx]
    pi.Y0Z1WWXX=pi.Y0Z1WXu[oo.wx]
    pi.Y1Z1WWXX=pi.Y1Z1WXu[oo.wx]



    XYZu<-expand.grid(Xu,c(0,1),c(0,1),Xu)
    XYZu=XYZu[,1]*100+XYZu[,2]*10+XYZu[,3]
    XXYYZZ<-c(X[delta==1],XX)*100+c(Y[delta==1],YY)*10+c(Z[delta==1],ZZ)
    oo.xyz<-match(XXYYZZ,XYZu)
    pi.XXYYZZ<-pi.XYZ[oo.xyz]


    #Xm<-as.matrix(cbind(1,c(Z[delta==1],ZZ),c(S1[delta==1],SS1),c(Z[delta==1],ZZ)*c(S1[delta==1],SS1),
    #c(W[delta==1],WW)-1,c(Z[delta==1],ZZ)*(c(W[delta==1],WW)-1)))

    Xm<-as.matrix(cbind(1,c(Z[delta==1],ZZ),c(S1[delta==1],SS1),c(Z[delta==1],ZZ)*c(S1[delta==1],SS1)))


    Y.long=as.vector(c(Y[delta==1],YY))

    iter<-0
    repeat{

    ### probability of delta=1 conditional on S1,W


    ### if random sample S from cases and controls in vaccine arm
    P.Z1<-pi.Y1Z1WWXX*R.probit(Z=1,S1=SS1,W=WW-1,beta=beta.old)+pi.Y0Z1WWXX*(1-R.probit(Z=1,S1=SS1,W=WW-1,beta=beta.old))


    P.Z0<-pi.Y0Z0WWXX*(1-R.probit(Z=0,S1=SS1,W=WW-1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ,S1=SS1,W=WW-1,beta=beta.old)
    ff<-ifelse(YY==1,ff,1-ff)

    ff<-ff/PP
    if (smooth==TRUE){
        ff<-ff*smoothWeights
    }


    if (arg==1) {
     if (length(ff)>0){
     denom<-tapply(ff,ind,sum)

     denom<-denom[ind]
     weights<-as.numeric(ff)/denom
     } else weights<-ff

    }


    fit=try.error(glm.fit(Xm,Y.long,weights=c(rep(1,sum(delta==1)),weights)/pi.XXYYZZ,family=binomial(link=probit),start=beta.old))

    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3


    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }
    #print(c(beta.new,iter))
    beta.old<-beta.new

    #print(iter)
    if (iter>500 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>500 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}



EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.X.Short<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,p.delta2,pd2.YZXW,beta){

  #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma



    Sval<-Sout$Sval
    lval2<-Sout$lval


    WX<-X*10+W
    WXu.ep<-expand.grid(Wu,Xu)
    WXu<-WXu.ep[,2]*10+WXu.ep[,1]

    #oo.2<-match(WX[delta1==1 & delta2==0],WXu)
    oo.2<-match(WX[delta1==1 & delta2==0], unlist(Sout$WXname))

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval2[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    XX.2<-rep(X[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(XX.2*10+WW.2,WXu)

    #############

    if (sum(delta1==0)>0){
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


    ind.1<-unlist(dataLong$indlong)
    YY.1<-unlist(dataLong$Ylong)
    WW.1<-unlist(dataLong$Wlong)
    ZZ.1<-unlist(dataLong$Zlong)
    XX.1<-unlist(dataLong$Xlong)
    SS1.1<-unlist(dataLong$Slong)
    Sweights.1<-unlist(dataLong$Sweightslong)

    oo.11<-match(XX.1*10+WW.1,WXu)
    } else {
      YY.1<-WW.1<-ZZ.1<-XX.1<-SS1.1<-Sweights.1<-oo.11<-integer(0)
    }

   # Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1),
#c(W[delta2==1],WW.2,WW.1)-1,c(Z[delta2==1],ZZ.2,ZZ.1)*(c(W[delta2==1],WW.2,WW.1)-1)
 #   ))

  Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))



    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2-1,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1-1,beta=beta.old))
    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))

    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    PP.delta2<-p.delta2[oo.11]

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1-1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff/PP*PP.delta2*Sweights.1



     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(Xm,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}





#########################################


EM.cc.CPV.Probit.Small.Nonpar2.SubW.X.Short<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,beta){

    Sval<-Sout$Sval
    lval2<-Sout$lval


    WX<-X*10+W
    WXu.ep<-expand.grid(Wu,Xu)
    WXu<-WXu.ep[,2]*10+WXu.ep[,1]

    #oo.2<-match(WX[delta1==1 & delta2==0],WXu)
    oo.2<-match(WX[delta1==1 & delta2==0], unlist(Sout$WXname))

    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval2[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    XX.2<-rep(X[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    SS1.2<-unlist(Sval[oo.2])
    Sweights.2<-unlist(Sout$Sweights[oo.2])
    oo.22<-match(XX.2*10+WW.2,WXu)

    #############

    if (sum(delta1==0)>0){
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


    ind.1<-unlist(dataLong$indlong)
    YY.1<-unlist(dataLong$Ylong)
    WW.1<-unlist(dataLong$Wlong)
    ZZ.1<-unlist(dataLong$Zlong)
    XX.1<-unlist(dataLong$Xlong)
    SS1.1<-unlist(dataLong$Slong)
    Sweights.1<-unlist(dataLong$Sweights12long)

    oo.11<-match(XX.1*10+WW.1,WXu)
    } else {
      YY.1<-WW.1<-ZZ.1<-XX.1<-SS1.1<-Sweights.1<-oo.11<-integer(0)
    }


  #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))

  #Xm<-as.matrix(cbind(1,c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1),
  #c(W[delta2==1],WW.2,WW.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(W[delta2==1],WW.2,WW.1)
  #))

  Xm<-as.matrix(cbind(1,c(Z[delta2==1],ZZ.2,ZZ.1), c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))



    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

   # P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))
#
#    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))
#    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff*Sweights.2
    #ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    #P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old))
#
#    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1-1,beta=beta.old))
#    PP<-P.Z1*mean(Z)+P.Z0*(1-mean(Z))
#
#    ### probability of delta2==1 conditional on delta1==1, W, and X
#
#    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
#    PP.delta2<-p.delta2[oo.11]

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff*Sweights.1
    #ff<-ff/PP*PP.delta2*Sweights.1



     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     denom<-denom[ind.1]
     weights.1<-as.numeric(ff)/denom
     } else weights.1<-ff


    fit=try.error(glm.fit(Xm,Y.long,weights=c(rep(1,sum(delta2==1)),weights.2,weights.1),family=binomial(link=probit),start=beta.old))


    out<-fit$coef

    beta.new<-out

    iter<-iter+1

    diffcount<-sum(abs(beta.new-beta.old))<1e-3

    #print(iter);
    #print(beta.new);
    if (is.na(diffcount)) break
    else{
     if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
    }

    beta.old<-beta.new

    #print(iter)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new))
}

######################

covDFcc.Nonpar2.X<-function(S, W, X, Wu, Xu,weights){

#########################################
### function to estimate covariate adjusted control distribution
## nonparametrically, conditional on W only
##### S: marker
##### Z: treatment
##### W, X: categorical baseline covariate
##### Wu, Xu: unique categories for W

nobs<-length(S)
SS<-ll<-Sweights<-WXname<-list()

index=0
for (j in 1:length(Xu)){
for (i in 1:length(Wu)){
     index=index+1
     oo<-order(S[W==Wu[i] & X==Xu[j]])
     SS[[index]]<-(S[W==Wu[i] & X==Xu[j]])[oo]
     ll[[index]]<-sum(W==Wu[i] & X==Xu[j])
     Sweights[[index]]<-(weights[W==Wu[i] & X==Xu[j]])[oo]
     WXname[[index]]<-Xu[j]*10+Wu[i]
}
}
return(list(Sval=SS,lval=ll,Sweights=Sweights,WXname=WXname))
}


covDFcc.Nonpar2.Full<-function(S, W, Wu){

#########################################
### function to estimate covariate adjusted control distribution
## nonparametrically, conditional on W only
##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
SS<-ll<-Sweights<-SSW<-list()
for (i in 1:length(Wu)){
 oo<-order(S[W==Wu[i]])
 SS[[i]]<-unique((S[W==Wu[i]])[oo])
 SSW[[i]]<-unique(SS[[i]]+Wu[i]*10)
 ll[[i]]<-length(unique(S[W==Wu[i]]))
 #Sweights[[i]]<-(weights[W==Wu[i]])[oo]
 Sweights[[i]]<-table(S[W==Wu[i]])/sum(W==Wu[i])
}
return(list(Sval=SS,SWval=SSW,lval=ll,Sweights=Sweights))
}


covDFcc.Nonpar2.Full.weights<-function(S, W, Wu,weights){

#########################################
### function to estimate covariate adjusted control distribution
## nonparametrically, conditional on W only
##### S: marker
##### Z: treatment
##### W: categorical baseline covariate
##### Wu: unique categories for W

nobs<-length(S)
SS<-ll<-Sweights<-SSW<-list()
for (i in 1:length(Wu)){
 oo<-order(S[W==Wu[i]])
 SS[[i]]<-unique((S[W==Wu[i]])[oo])
 SSW[[i]]<-unique(SS[[i]]+Wu[i]*10)
 ll[[i]]<-length(unique(S[W==Wu[i]]))
 #Sweights[[i]]<-(weights[W==Wu[i]])[oo]
 Sweights[[i]]<-tapply(weights[W==Wu[i]],S[W==Wu[i]],sum)/sum(weights[W==Wu[i]])
}
return(list(Sval=SS,SWval=SSW,lval=ll,Sweights=Sweights))
}




getDataLong<-function(Y,Z,X,Wu,Xu,Sout,Wout){


Slong<-Sweights1long<-Sweights12long<-Wlong<-Ylong<-Zlong<-Xlong<-indlong<-list()


for (i in 1:length(X)){
  ooX=match(X[i],Xu)
  Wval=unlist(Wout$Sval[ooX])
  WXind=match(X[i]*10+Wval,unlist(Sout$WXname))
  Sweights1long[[i]]=rep(unlist(Wout$Sweights[ooX]),unlist(Sout$lval[WXind]))
  Sweights12long[[i]]=Sweights1long[[i]]*(unlist(Sout$Sweights[WXind]))
  Wlong[[i]]=rep(Wval,unlist(Sout$lval[WXind]))
  Slong[[i]]<-Sout$Sval[WXind]

  indlong[[i]]<-rep(i,sum(unlist(Sout$lval[WXind])))
  Ylong[[i]]<-rep(Y[i],sum(unlist(Sout$lval[WXind])))
  Zlong[[i]]<-rep(Z[i],sum(unlist(Sout$lval[WXind])))
  Xlong[[i]]<-rep(X[i],sum(unlist(Sout$lval[WXind])))
}
return(list(Wlong=Wlong,Ylong=Ylong,Xlong=Xlong,Zlong=Zlong,Slong=Slong,indlong=indlong,
Sweights1long=Sweights1long,Sweights12long=Sweights12long))
}


###################


EM.cc.CPV.Probit.Small.Nonpar2.weighted<-function(Z,Sout,S1,W,Y,Wu,delta,beta,weights1){


Sval<-Sout$Sval
lval<-Sout$lval
Sweights<-Sout$Sweights


oo<-match(W[delta==0],Wu)
ind<-rep(1:sum(delta==0), unlist(lval[oo]))

YY<-rep(Y[delta==0],unlist(lval[oo]))

WW<-rep(W[delta==0],unlist(lval[oo]))

WWeights1<-rep(weights1[delta==0],unlist(lval[oo]))

ZZ<-rep(Z[delta==0],unlist(lval[oo]))


SS1<-unlist(Sval[oo])

resWeights<-unlist(Sweights[oo])

WWeights1.long<-c(weights1[delta==1],WWeights1)


beta.old<-beta


iter<-0
repeat{

ff<-R.probit(Z=ZZ,S1=SS1,W=WW,beta=beta.old)
ff<-ifelse(YY==1,ff,1-ff)
ff<-ff*resWeights


 if (length(ff)>0){
 denom<-tapply(ff,ind,sum)
 denom<-denom[ind]
 weights<-as.numeric(ff)/denom
 } else weights<-ff

dataout<-data.frame(Y=c(Y[delta==1],YY),Z=c(Z[delta==1],ZZ),W=c(W[delta==1],WW),
S1=c(S1[delta==1],SS1),weights=c(rep(1,sum(delta==1)),weights))

fit<-glm(Y~Z+S1+I(Z*S1),weights=weights/WWeights1.long,data=dataout,family=binomial(link=probit),start=beta)


out<-fit$coef

beta.new<-out

iter<-iter+1

diffcount<-sum(abs(beta.new-beta.old))<1e-3

if (is.na(diffcount)) break
else{
 if (sum(abs(beta.new-beta.old))<1e-3 & sum(abs(beta.new-beta.old))<sum(abs(beta.old))*1e-3) break
}

beta.old<-beta.new


if (iter>200 | is.na(diffcount)) break
}
loglike<-logLik(fit)
if (iter>200 | is.na(diffcount)) {
  beta.new<-rep(NA,length(beta.old))
  loglike<-NA
}
#return(c(beta.new,iter,loglike))
return(c(beta.new))
}
