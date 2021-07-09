
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


covDFcc.Nonpar2.X.New<-function(S, W, X, Wu, Xu,weights){

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
     Sweights[[index]]<-(weights[W==Wu[i] & X==Xu[j]])[oo]/sum((weights[W==Wu[i] & X==Xu[j]])[oo])
     WXname[[index]]<-Xu[j]*10+Wu[i]
}
}
return(list(Sval=SS,lval=ll,Sweights=Sweights,WXname=WXname))
}




getDummy.COVID.M<-function(dat,arg=0){
    
  dat$MinorityInd=floor(dat$W/100)==1
  dat$HighriskInd=floor((dat$W-dat$MinorityInd*100)/10)==1
  dat$Senior=dat$W-dat$MinorityInd*100-dat$HighriskInd*10==1

  return(dat)
}



getDummy.COVID.E<-function(dat){
### using baseline covariates as BIP
    
  dat$RM1=floor(dat$W/100)==1
  dat$RM2=floor(dat$W/100)==2
  dat$RM3=floor(dat$W/100)==3
  dat$HighriskInd=floor((dat$W-floor(dat$W/100)*100)/10)==1
  dat$Senior=dat$W-floor(dat$W/100)*100-dat$HighriskInd*10
  
  return(dat)
}

getDummy.COVID.EB<-function(dat){
### assuming BIP in case-cohort samples
    
  dat$RM1=floor(dat$X/100)==1
  dat$RM2=floor(dat$X/100)==2
  dat$RM3=floor(dat$X/100)==3
  dat$HighriskInd=floor((dat$X-floor(dat$X/100)*100)/10)==1
  dat$Senior=dat$X-floor(dat$X/100)*100-dat$HighriskInd*10
  
  
  dat$W2<-dat$W==2
  dat$W3<-dat$W==3
  dat$W4<-dat$W==4
  return(dat)
}



getDummy.COVID.EB.WX<-function(dat){
   
  dat$X=floor(dat$WX/10)
  dat$W=dat$WX-dat$X*10  
  dat$RM1=floor(dat$X/100)==1
  dat$RM2=floor(dat$X/100)==2
  dat$RM3=floor(dat$X/100)==3
  dat$HighriskInd=floor((dat$X-floor(dat$X/100)*100)/10)==1
  dat$Senior=dat$X-floor(dat$X/100)*100-dat$HighriskInd*10
  
  
  dat$W2<-dat$W==2
  dat$W3<-dat$W==3
  dat$W4<-dat$W==4
  return(dat)
}

R.probit.COVID.M<-function(Z,S1,W,X,varlist,beta){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.COVID.M(datin,arg=arg)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}

R.probit.COVID.E<-function(Z,S1,W,X,varlist,beta){

### Use baseline covariates as BIP
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.COVID.E(datin)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}


R.probit.COVID.EB<-function(Z,S1,W,X,varlist,beta){

### assume BIP in case-cohort sample
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.COVID.EB(datin)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}

R.probit.COVID.EB.WX<-function(Z,S1,WX,varlist,beta){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,WX=WX)
   datin=getDummy.COVID.EB.WX(datin)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}


#### X is 1, W includes a set of baseline covariates and stratum measured from all trial participants

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.M<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta){

   #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma
 
    PZ=mean(Z)

    Sval<-Sout$Sval
    lval2<-Sout$lval


    WX<-X*10+W
    WXu.ep<-expand.grid(Wu,Xu)
    WXu<-WXu.ep[,2]*10+WXu.ep[,1]

    ### generate augmented dataset for those with delta1==1 & delta2==0
    #oo.2<-match(WX[delta1==1 & delta2==0],WXu)
    oo.2<-match(WX[delta1==1 & delta2==0], unlist(Sout$WXname))

    # (WX[delta1==1 & delta2==0])[1:5]
    #     [1] 43 11 51 41 42
    #    > unlist((Sout$WXname))[oo.2[1:5]]
    #    [1] 43 11 51 41 42


    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval2[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    XX.2<-rep(X[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(XX.2*10+WW.2,WXu)

    #############

    if (sum(delta1==0)>0){
    ### generate augmented data for those with delta1=0
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


#> getDataLong
#function(Y,Z,X,Wu,Xu,Sout,Wout){
#
#Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
#
#Slong<-Sweightslong<-Wlong<-Ylong<-Zlong<-Xlong<-indlong<-list()
#
#for (i in 1:length(X)){
#  ooX=match(X[i],Xu)
#  Wval=unlist(Wout$Sval[ooX])
#  WXind=match(X[i]*10+Wval,unlist(Sout$WXname))
#  Sweightslong[[i]]=rep(unlist(Wout$Sweights[ooX]),unlist(Sout$lval[WXind]))
#  Wlong[[i]]=rep(Wval,unlist(Sout$lval[WXind]))
#  Slong[[i]]<-Sout$Sval[WXind]
#
#  indlong[[i]]<-rep(i,sum(unlist(Sout$lval[WXind])))
#  Ylong[[i]]<-rep(Y[i],sum(unlist(Sout$lval[WXind])))
#  Zlong[[i]]<-rep(Z[i],sum(unlist(Sout$lval[WXind])))
#  Xlong[[i]]<-rep(X[i],sum(unlist(Sout$lval[WXind])))
#}
#return(list(Wlong=Wlong,Ylong=Ylong,Xlong=Xlong,Zlong=Zlong,Slong=Slong,indlong=indlong,Sweightslong=Sweightslong))
#}



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

   # Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1),
#c(W[delta2==1],WW.2,WW.1)-1,c(Z[delta2==1],ZZ.2,ZZ.1)*(c(W[delta2==1],WW.2,WW.1)-1)
 #   ))

    dat.Xm<-data.frame(Z=c(Z[delta2==1],ZZ.2,ZZ.1),S1=c(S1[delta2==1],SS1.2,SS1.1),W=c(W[delta2==1],WW.2,WW.1),X=c(X[delta2==1],XX.2,XX.1))
    dat.Xm<-getDummy.COVID.M(dat.Xm)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.COVID.M(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit.COVID.M(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.COVID.M(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.COVID.M(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     #denom<-denom[ind.2] ## only works if there is no gap in ind.2
     denom<-denom[match(ind.2,names(denom))]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.COVID.M(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit.COVID.M(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.COVID.M(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.COVID.M(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    #ff<-ff/PP*PP.delta2*Sweights.1
    ff<-ff/PP*Sweights.1


     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     #denom<-denom[ind.1]  ## only works if there is no gap in ind.1
     denom<-denom[match(ind.1,names(denom))]
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
    #print(beta.new)
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





#### X is 1, W includes a set of baseline covariates and stratum measured from all trial participants

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.E<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta){

   ### assume BIP available from case-cohort sample
   #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma
 
    PZ=mean(Z)

    Sval<-Sout$Sval
    lval2<-Sout$lval


    WX<-X*10+W
    WXu.ep<-expand.grid(Wu,Xu)
    WXu<-WXu.ep[,2]*10+WXu.ep[,1]

    ### generate augmented dataset for those with delta1==1 & delta2==0
    #oo.2<-match(WX[delta1==1 & delta2==0],WXu)
    oo.2<-match(WX[delta1==1 & delta2==0], unlist(Sout$WXname))

    # (WX[delta1==1 & delta2==0])[1:5]
    #     [1] 43 11 51 41 42
    #    > unlist((Sout$WXname))[oo.2[1:5]]
    #    [1] 43 11 51 41 42


    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval2[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    XX.2<-rep(X[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(XX.2*10+WW.2,WXu)

    #############

    if (sum(delta1==0)>0){
    ### generate augmented data for those with delta1=0
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


#> getDataLong
#function(Y,Z,X,Wu,Xu,Sout,Wout){
#
#Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
#
#Slong<-Sweightslong<-Wlong<-Ylong<-Zlong<-Xlong<-indlong<-list()
#
#for (i in 1:length(X)){
#  ooX=match(X[i],Xu)
#  Wval=unlist(Wout$Sval[ooX])
#  WXind=match(X[i]*10+Wval,unlist(Sout$WXname))
#  Sweightslong[[i]]=rep(unlist(Wout$Sweights[ooX]),unlist(Sout$lval[WXind]))
#  Wlong[[i]]=rep(Wval,unlist(Sout$lval[WXind]))
#  Slong[[i]]<-Sout$Sval[WXind]
#
#  indlong[[i]]<-rep(i,sum(unlist(Sout$lval[WXind])))
#  Ylong[[i]]<-rep(Y[i],sum(unlist(Sout$lval[WXind])))
#  Zlong[[i]]<-rep(Z[i],sum(unlist(Sout$lval[WXind])))
#  Xlong[[i]]<-rep(X[i],sum(unlist(Sout$lval[WXind])))
#}
#return(list(Wlong=Wlong,Ylong=Ylong,Xlong=Xlong,Zlong=Zlong,Slong=Slong,indlong=indlong,Sweightslong=Sweightslong))
#}

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

   # Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1),
#c(W[delta2==1],WW.2,WW.1)-1,c(Z[delta2==1],ZZ.2,ZZ.1)*(c(W[delta2==1],WW.2,WW.1)-1)
 #   ))

    dat.Xm<-data.frame(Z=c(Z[delta2==1],ZZ.2,ZZ.1),S1=c(S1[delta2==1],SS1.2,SS1.1),W=c(W[delta2==1],WW.2,WW.1),X=c(X[delta2==1],XX.2,XX.1))
    dat.Xm<-getDummy.COVID.E(dat.Xm)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.COVID.E(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit.COVID.E(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.COVID.E(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.COVID.E(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP

     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     #denom<-denom[ind.2] ## only works if there is no gap in ind.2
     denom<-denom[match(ind.2,names(denom))]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.COVID.E(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit.COVID.E(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.COVID.E(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]

    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.COVID.E(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    #ff<-ff/PP*PP.delta2*Sweights.1
    ff<-ff/PP*Sweights.1

     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     #denom<-denom[ind.1]  ## only works if there is no gap in ind.1
     denom<-denom[match(ind.1,names(denom))]
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
    #print(beta.new)
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




#### X includes confounders and sampling strata, W (categoried as W2, W3, W4) measured in a subset of participants with X measured

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID.EB<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta){

   ### assume BIP available from case-cohort sample
   #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma
 
    PZ=mean(Z)

    Sval<-Sout$Sval
    lval2<-Sout$lval


    WX<-X*10+W
    WXu.ep<-expand.grid(Wu,Xu)
    WXu<-WXu.ep[,2]*10+WXu.ep[,1]

    ### generate augmented dataset for those with delta1==1 & delta2==0
    #oo.2<-match(WX[delta1==1 & delta2==0],WXu)
    oo.2<-match(WX[delta1==1 & delta2==0], unlist(Sout$WXname))

    # (WX[delta1==1 & delta2==0])[1:5]
    #     [1] 43 11 51 41 42
    #    > unlist((Sout$WXname))[oo.2[1:5]]
    #    [1] 43 11 51 41 42


    ind.2<-rep(1:sum(delta1==1 & delta2==0), unlist(lval2[oo.2]))
    YY.2<-rep(Y[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    WW.2<-rep(W[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    ZZ.2<-rep(Z[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    XX.2<-rep(X[delta1==1 & delta2==0],unlist(lval2[oo.2]))
    SS1.2<-unlist(Sval[oo.2])

    oo.22<-match(XX.2*10+WW.2,WXu)

    #############

    if (sum(delta1==0)>0){
    ### generate augmented data for those with delta1=0
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


#> getDataLong
#function(Y,Z,X,Wu,Xu,Sout,Wout){
#
#Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
#
#Slong<-Sweightslong<-Wlong<-Ylong<-Zlong<-Xlong<-indlong<-list()
#
#for (i in 1:length(X)){
#  ooX=match(X[i],Xu)
#  Wval=unlist(Wout$Sval[ooX])
#  WXind=match(X[i]*10+Wval,unlist(Sout$WXname))
#  Sweightslong[[i]]=rep(unlist(Wout$Sweights[ooX]),unlist(Sout$lval[WXind]))
#  Wlong[[i]]=rep(Wval,unlist(Sout$lval[WXind]))
#  Slong[[i]]<-Sout$Sval[WXind]
#
#  indlong[[i]]<-rep(i,sum(unlist(Sout$lval[WXind])))
#  Ylong[[i]]<-rep(Y[i],sum(unlist(Sout$lval[WXind])))
#  Zlong[[i]]<-rep(Z[i],sum(unlist(Sout$lval[WXind])))
#  Xlong[[i]]<-rep(X[i],sum(unlist(Sout$lval[WXind])))
#}
#return(list(Wlong=Wlong,Ylong=Ylong,Xlong=Xlong,Zlong=Zlong,Slong=Slong,indlong=indlong,Sweightslong=Sweightslong))
#}

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

   # Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1),
#c(W[delta2==1],WW.2,WW.1)-1,c(Z[delta2==1],ZZ.2,ZZ.1)*(c(W[delta2==1],WW.2,WW.1)-1)
 #   ))

    dat.Xm<-data.frame(Z=c(Z[delta2==1],ZZ.2,ZZ.1),S1=c(S1[delta2==1],SS1.2,SS1.1),W=c(W[delta2==1],WW.2,WW.1),X=c(X[delta2==1],XX.2,XX.1))
    dat.Xm<-getDummy.COVID.E(dat.Xm)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.COVID.EB(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit.COVID.EB(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.COVID.EB(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.COVID.EB(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP

     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     #denom<-denom[ind.2] ## only works if there is no gap in ind.2
     denom<-denom[match(ind.2,names(denom))]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.COVID.EB(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit.COVID.EB(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.COVID.EB(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]

    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.COVID.EB(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    #ff<-ff/PP*PP.delta2*Sweights.1
    ff<-ff/PP*Sweights.1

     if (length(ff)>0){
     denom<-tapply(ff,ind.1,sum)
     #denom<-denom[ind.1]  ## only works if there is no gap in ind.1
     denom<-denom[match(ind.1,names(denom))]
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
    #print(beta.new)
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



integVE.COVID.M<-function(Su,datin,varlist,beta){
   
   datin=datin[!is.na(datin$S1),]  ### Need this since might end up with W categories with no non-missing S(1)
   Wu=unique(sort(datin$W))
   datnew<-data.frame(S1=Su)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(.25,.5,.75),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(1/3,2/3),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(.5),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))

    fx.s=predict(fit,datnew,type='probs')
     PY1.SZ0<-PY1.SZ1<-PY1.SZ0.b<-PY1.SZ1.b<-rep(0,length(Su))

    for (i in 1:length(Wu)){

      ### compute P(Y=1,Z=1|S) and P(Y=1,Z=0|S) for each Su value in the population#####################

      PY1.SZ1=PY1.SZ1+R.probit.COVID.M(Z=1,S1=Su,W=Wu[i],X=1,varlist=varlist,beta)*fx.s[,i]
      PY1.SZ0=PY1.SZ0+R.probit.COVID.M(Z=0,S1=Su,W=Wu[i],X=1,varlist=varlist,beta)*fx.s[,i]
      

    }

###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
   return(list(Su=Su,VE=VE.S))
   }



integVE.COVID.E<-function(Su,datin,varlist,beta){
   ### use baseline covariates as BIP
 
  datin=datin[!is.na(datin$S1),]  ### Need this since might end up with W categories with no non-missing S(1)
   Wu=unique(sort(datin$W))
   datnew<-data.frame(S1=Su)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(.25,.5,.75),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(1/3,2/3),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(W~ns(S1,knots=quantile(S1,c(.5),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))

    fx.s=predict(fit,datnew,type='probs')
     PY1.SZ0<-PY1.SZ1<-PY1.SZ0.b<-PY1.SZ1.b<-rep(0,length(Su))

    for (i in 1:length(Wu)){

      ### compute P(Y=1,Z=1|S) and P(Y=1,Z=0|S) for each Su value in the population#####################

      PY1.SZ1=PY1.SZ1+R.probit.COVID.E(Z=1,S1=Su,W=Wu[i],X=1,varlist=varlist,beta)*fx.s[,i]
      PY1.SZ0=PY1.SZ0+R.probit.COVID.E(Z=0,S1=Su,W=Wu[i],X=1,varlist=varlist,beta)*fx.s[,i]      

    }

###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
   return(list(Su=Su,VE=VE.S))

}


integVE.COVID.EB<-function(Su,datin,varlist,beta){
   ### assume BIP available in case-cohort sample
   datin$WX<-datin$X*10+datin$W

   datin=datin[!is.na(datin$S1),]  ### Need this since might end up with W categories with no non-missing S(1)
   WXu=unique(sort(datin$WX))
   datnew<-data.frame(S1=Su)
    fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(.25,.5,.75),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(1/3,2/3),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))
    if (fit$convergence==1)
    fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(.5),na.rm=T),Boundary.knots=range(S1,na.rm=T)),weights=weights,data=datin))

    fx.s=predict(fit,datnew,type='probs')
     PY1.SZ0<-PY1.SZ1<-PY1.SZ0.b<-PY1.SZ1.b<-rep(0,length(Su))

    for (i in 1:length(WXu)){

      ### compute P(Y=1,Z=1|S) and P(Y=1,Z=0|S) for each Su value in the population#####################

      PY1.SZ1=PY1.SZ1+R.probit.COVID.EB.WX(Z=1,S1=Su,WX=WXu[i],varlist=varlist,beta)*fx.s[,i]
      PY1.SZ0=PY1.SZ0+R.probit.COVID.EB.WX(Z=0,S1=Su,WX=WXu[i],varlist=varlist,beta)*fx.s[,i]
      
    }

###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
   return(list(Su=Su,VE=VE.S))
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
