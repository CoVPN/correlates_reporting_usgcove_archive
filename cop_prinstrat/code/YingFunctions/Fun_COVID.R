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

## old version
#getDummy.New<-function(dat,arg=0){
#  if (arg==0) {
#
#  dat$VB=floor(dat$W/100)==1
#  dat$RS=floor((dat$W-floor(dat$W/100)*100)/10)==1
#  
#  dat$W2=dat$W-floor(dat$W/100)*100-floor(dat$W/10)*10==2
#  dat$W3=dat$W-floor(dat$W/100)*100-floor(dat$W/10)*10==3
#  dat$W4=dat$W-floor(dat$W/100)*100-floor(dat$W/10)*10==4
#  } else {
#    WX=dat$X*10+dat$W
#    dat$VB=floor(WX/100)==1
#    dat$RS=floor((WX-floor(WX/100)*100)/10)==1
#  
#    dat$W2=WX-floor(WX/100)*100-floor(WX/10)*10==2
#    dat$W3=WX-floor(WX/100)*100-floor(WX/10)*10==3
#    dat$W4=WX-floor(WX/100)*100-floor(WX/10)*10==4
#  
#  }
#  return(dat)
#}

### Modified 06/01/2021
getDummy.New<-function(dat,arg=0){
  if (arg==0) {

  dat$VB=floor(dat$W/100)==1
  dat$RS=floor((dat$W-floor(dat$W/100)*100)/10)==1
  
  dat$W2=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==2
  dat$W3=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==3
  dat$W4=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==4
  } else {
    WX=dat$X*10+dat$W
    dat$VB=floor(WX/100)==1
    dat$RS=floor((WX-floor(WX/100)*100)/10)==1
  
    dat$W2=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==2
    dat$W3=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==3
    dat$W4=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==4
  
  }
  return(dat)
}
getDummy.New.X<-function(dat){

  dat$VB=floor(dat$X/1000)==1
  dat$RS=floor((dat$X-floor(dat$X/1000)*1000)/100)==1
  
  dat$W2=dat$W==2
  dat$W3=dat$W==3
  dat$W4=dat$W==4
  
  return(dat)
}
getDummy.COVID<-function(dat,arg=0){
  
  dat$MinorityInd=floor(dat$W/100)==1
  dat$HighriskInd=floor((dat$W-floor(dat$W/100)*100)/10)==1
  dat$Agegr=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==1

  return(dat)
}


getDummy.COVID<-function(dat,arg=0){
  
  if (arg==0){
  dat$MinorityInd=floor(dat$W/100)==1
  dat$HighriskInd=floor((dat$W-dat$Minority*100)/10)==1
  dat$Agegr=dat$W-dat$Minority*100-dat$HighriskInd*10==1
  } {
  
    WX=dat$X*10+dat$W
    dat$MinorityInd=floor(WX/1000)==1
    dat$HighriskInd=floor((WX-dat$MinorityInd*1000)/100)==1
    dat$Agegr=floor((WX-dat$MinorityInd*1000-dat$HighriskInd*100)/10)==1
    
    
    dat$W2=WX-dat$MinorityInd*1000-dat$HighriskInd*100-dat$Agegr*10==2
    dat$W3=WX-dat$MinorityInd*1000-dat$HighriskInd*100-dat$Agegr*10==3
    dat$W4=WX-dat$MinorityInd*1000-dat$HighriskInd*100-dat$Agegr*10==4
  }
  return(dat)
}

### Modified 06/01/2021
getDummy.New<-function(dat,arg=0){
  if (arg==0) {

  dat$VB=floor(dat$W/100)==1
  dat$RS=floor((dat$W-floor(dat$W/100)*100)/10)==1
  
  dat$W2=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==2
  dat$W3=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==3
  dat$W4=dat$W-floor(dat$W/100)*100-floor((dat$W-floor(dat$W/100)*100)/10)*10==4
  } else {
    WX=dat$X*10+dat$W
    dat$VB=floor(WX/100)==1
    dat$RS=floor((WX-floor(WX/100)*100)/10)==1
  
    dat$W2=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==2
    dat$W3=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==3
    dat$W4=WX-floor(WX/100)*100-floor((WX-floor(WX/100)*100)/10)*10==4
  
  }
  return(dat)
}

R.probit.in<-function(Z,S1,W,X,varlist,beta,arg=0){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.New(datin,arg=arg)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
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

R.probit.COVID.E.WX<-function(Z,S1,WX,varlist,beta){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,WX=WX)
   datin=getDummy.COVID.E.WX(datin)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}

R.probit.COVID<-function(Z,S1,W,X,varlist,beta,arg=0){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.COVID(datin,arg=arg)
   oo<-match(varlist,names(datin))
   X<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(X%*%as.matrix(beta))
   }
   return(Risk)
}

R.probit.in.X<-function(Z,S1,W,X,varlist,beta){
####### calculate the logit risk ############
   Risk<-NULL
   if (length(S1)>0){
   datin=data.frame(Z=Z,S1=S1,W=W,X=X)
   datin=getDummy.New.X(datin)
   oo<-match(varlist,names(datin))
   XX<-as.matrix(cbind(1,datin$Z,datin$S1,datin$S1*datin$Z,datin[,oo]))
   Risk<-pnorm(XX%*%as.matrix(beta))
   }
   return(Risk)
} 

#Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;beta=c(0,0,0,0);varlist=NULL;arg=0
EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta,arg=0){

   #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma
   #Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;varlist=NULL;beta=c(0,0,0,0)
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
    dat.Xm<-getDummy.New(dat.Xm,arg=arg)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.in(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg)+
        pd2.YZXW[1,2,oo.22]*(1-R.probit.in(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.in(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.in(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.in(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg)+
        pd2.YZXW[1,2,oo.11]*(1-R.probit.in(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.in(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.in(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg)
    ff<-ifelse(YY.1==1,ff,1-ff)
    #ff<-ff/PP*PP.delta2*Sweights.1
    ff<-ff/PP*Sweights.1


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


  integVE<-function(Su,Wu,datin,varlist,beta,arg=0){
   
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

      PY1.SZ1=PY1.SZ1+R.probit.in(Z=1,S1=Su,W=Wu[i],X=1,varlist=varlist,beta,arg=arg)*fx.s[,i]
      PY1.SZ0=PY1.SZ0+R.probit.in(Z=0,S1=Su,W=Wu[i],X=1,varlist=varlist,beta,arg=arg)*fx.s[,i]
      

    }

###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
   return(list(Su=Su,VE=VE.S))
   }


#Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;beta=c(0,0,0,0);varlist=NULL

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.2020.X<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta){

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
    dat.Xm<-getDummy.New.X(dat.Xm)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.in.X(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit.in.X(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.in.X(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.in.X(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     #denom<-denom[ind.2] ## only works if there is no gap in ind.2
     denom<-denom[match(ind.2,names(denom))]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.in.X(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit.in.X(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.in.X(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.in.X(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old)
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

  
 integVE.X<-function(Su,varlist=NULL,datin,beta){

   datnew<-data.frame(S1=Su)
    getDummy.X(datin)
    if (!is.null(varlist)){
        fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(.25,.5,.75)),Boundary.knots=range(S1)),weights=weights,data=datin))
        if (fit$convergence==1)
        fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(1/3,2/3)),Boundary.knots=range(S1)),weights=weights,data=datin))
        if (fit$convergence==1)
        fit=try.error(multinom(WX~ns(S1,knots=quantile(S1,c(.5)),Boundary.knots=range(S1)),weights=weights,data=datin))
    
        fx.s=predict(fit,datnew,type='probs')
         PY1.SZ0<-PY1.SZ1<-PY1.SZ0.b<-PY1.SZ1.b<-rep(0,length(Su))
    
        for (i in 1:length(WXu)){
          ### compute P(Y=1,Z=1|S) and P(Y=1,Z=0|S) for each Su value in the population#####################
    
          PY1.SZ1=PY1.SZ1+R.probit.WX(Z=1,S1=Su,W=WXu[i],beta)*fx.s[,i]
          PY1.SZ0=PY1.SZ0+R.probit.WX(Z=0,S1=Su,W=WXu[i],beta)*fx.s[,i]
          
        }

    } else {
      PY1.SZ1=R.probit(Z=1,S1=Su,W=NULL,beta)
      PY1.SZ0=R.probit(Z=0,S1=Su,W=NULL,beta)
    } 
    
###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
    
   return(list(Su=Su,VE=VE.S))
   }


EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New.COVID<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,varlist=NULL,beta,arg=0){

   #Sout=Sout.NC.All.X;Wout=Wout.NC.All.X;beta=gamma
   #Sout=Sout.NC.All.X;Wout=WoutA.NC.All.X;varlist=NULL;beta=c(0,0,0,0)
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
    dat.Xm<-getDummy.COVID(dat.Xm,arg=arg)
    #Xm<-as.matrix(cbind(1,Z=c(Z[delta2==1],ZZ.2,ZZ.1),c(S1[delta2==1],SS1.2,SS1.1),c(Z[delta2==1],ZZ.2,ZZ.1)*c(S1[delta2==1],SS1.2,SS1.1)))
    Xm<-as.matrix(cbind(1,dat.Xm$Z,dat.Xm$S1,dat.Xm$Z*dat.Xm$S1,dat.Xm[,match(varlist,names(dat.Xm))]))

    Y.long=as.vector(c(Y[delta2==1],YY.2,YY.1))


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit.COVID(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg)+
        pd2.YZXW[1,2,oo.22]*(1-R.probit.COVID(Z=1,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit.COVID(Z=0,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

    ff<-R.probit.COVID(Z=ZZ.2,S1=SS1.2,W=WW.2,X=XX.2,varlist=varlist,beta=beta.old,arg=arg)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit.COVID(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg)+
        pd2.YZXW[1,2,oo.11]*(1-R.probit.COVID(Z=1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg))

    P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit.COVID(Z=0,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg))
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit.COVID(Z=ZZ.1,S1=SS1.1,W=WW.1,X=XX.1,varlist=varlist,beta=beta.old,arg=arg)
    ff<-ifelse(YY.1==1,ff,1-ff)
    #ff<-ff/PP*PP.delta2*Sweights.1
    ff<-ff/PP*Sweights.1


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


   integVE.COVID<-function(Su,Wu,datin,varlist,beta,arg=0){
   
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

      PY1.SZ1=PY1.SZ1+R.probit.COVID(Z=1,S1=Su,W=Wu[i],X=1,varlist=varlist,beta,arg=arg)*fx.s[,i]
      PY1.SZ0=PY1.SZ0+R.probit.COVID(Z=0,S1=Su,W=Wu[i],X=1,varlist=varlist,beta,arg=arg)*fx.s[,i]
      

    }

###################### VE(S) curve ##################
    VE.S<-1-PY1.SZ1/PY1.SZ0  ## in the population
   return(list(Su=Su,VE=VE.S))
   }
