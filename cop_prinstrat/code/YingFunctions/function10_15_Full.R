
EM.cc.CPV.Probit.Small.Nonpar2.NC.Weighted.X.Short.New<-function(Z,Sout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd1.YZX,pd2.YZXW,beta){

   #Sout=Sout.NC.All.X;pd1.YZX=pd1.YZX.true;pd2.YZXW=pd2.YZXW.true;beta=gamma;
 

      PZ=mean(Z)

        ## probability of delta_1=1 given Y=y,Z=z,X=x
       # pi.Y0Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,Xu)
#        pi.Y1Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,Xu)
#        pi.Y0Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,Xu)
#        pi.Y1Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,Xu)


        pi.Y0Z0Xu<-pd1.YZX[1,1,]
        pi.Y1Z0Xu<-pd1.YZX[2,1,]
        pi.Y0Z1Xu<-pd1.YZX[1,2,]
        pi.Y1Z1Xu<-pd1.YZX[2,2,]

        pi.XYZ<-c(pi.Y0Z0Xu,pi.Y1Z0Xu,pi.Y0Z1Xu,pi.Y1Z1Xu)


        WXu.ep<-expand.grid(Wu,Xu)


       #    pi.Y1Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,1)
#           pi.Y0Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,1)
#           pi.Y1Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,0)
#           pi.Y0Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,0)

       pi.Y1Z1WXu<-pd2.YZXW[2,2,]
       pi.Y0Z1WXu<-pd2.YZXW[1,2,]
       pi.Y1Z0WXu<-pd2.YZXW[2,1,]
       pi.Y0Z0WXu<-pd2.YZXW[1,1,]



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



    oo<-match(WX[delta==0],WXu)

    ind<-rep(1:sum(delta==0), unlist(lval[oo]))


    YY<-rep(Y[delta==0],unlist(lval[oo]))


    WW<-rep(W[delta==0],unlist(lval[oo]))


    ZZ<-rep(Z[delta==0],unlist(lval[oo]))

    XX<-rep(X[delta==0],unlist(lval[oo]))


    SS1<-unlist(Sval[oo])

    WWXX<-WW+XX*10


   #### probability of P(Z=z|X=x)
   #PZ.X=table(X,Z)/as.numeric(table(X))

   oo.x=match(XX,Xu)

    beta.old<-beta



    oo.wx<-match(WWXX,WXu)
    pi.Y0Z0WWXX=pi.Y0Z0WXu[oo.wx]
    pi.Y1Z0WWXX=pi.Y1Z0WXu[oo.wx]
    pi.Y0Z1WWXX=pi.Y0Z1WXu[oo.wx]
    pi.Y1Z1WWXX=pi.Y1Z1WXu[oo.wx]



    XYZu<-expand.grid(Xu,c(0,1),c(0,1))
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
    PP<-P.Z1*PZ+P.Z0*(1-PZ)


    #PP=P.Z1*PZ.X[oo.x,2]+P.Z0*PZ.X[oo.x,1]

    ff<-R.probit(Z=ZZ,S1=SS1,W=WW-1,beta=beta.old)
    ff<-ifelse(YY==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind,sum)

     denom<-denom[ind]
     weights<-as.numeric(ff)/denom
     } else weights<-ff

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
    print(beta.new)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new,iter))
}



####################3

EM.cc.CPV.Probit.Small.Nonpar2.NC.SubW.X.Short.New<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,p.delta2,pd2.YZXW,beta){

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
##Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
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


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

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
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


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


EM.cc.CPV.Probit.Small.Nonpar2.NC.SubWA.X.Short.New<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,beta){

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
##Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
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


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x2=match(XX.2,Xu)
    #PP=P.Z1*PZ.X[oo.x2,2]+P.Z0*PZ.X[oo.x2,1]

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
    PP<-P.Z1*PZ+P.Z0*(1-PZ)

    #oo.x1=match(XX.1,Xu)
    #PP=P.Z1*PZ.X[oo.x1,2]+P.Z0*PZ.X[oo.x1,1]


    ### probability of delta2==1 conditional on delta1==1, W, and X

    #PP.delta2<-rep(p.delta2,sum(delta1==0)*unlist(lval))
    #PP.delta2<-p.delta2[oo.11]

    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1-1,beta=beta.old)
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
    print(beta.new)
    if (iter>200 | is.na(diffcount)) break
    }
    #loglike<-logLik(fit)
    loglike<-fit$deviance
    if (iter>200 | is.na(diffcount)) {
      beta.new<-rep(NA,length(beta.old))
      loglike<-NA
    }
    #return(c(beta.new,iter,loglike))
    return(c(beta.new,iter))
}




############# Estimated likelihood ###################




EM.cc.CPV.Probit.Small.Nonpar2.EL.SubW.X.Short.New<-function(Z,Sout,Wout,S1,W,Y,X,Wu,Xu,delta1,delta2,beta){

   #Sout=Sout.EL.All.X;Wout=Wout.EL.All.X;beta=gamma

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

    Sweights.2<-unlist(Sout$Sweights[oo.2])

    oo.22<-match(XX.2*10+WW.2,WXu)



    #############

    if (sum(delta1==0)>0){
    ### generate augmented data for those with delta1=0
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout,Wout)


#> getDataLong
#function(Y,Z,X,Wu,Xu,Sout,Wout){
#
##Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
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


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    #P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))

    #P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))


    #PP<-P.Z1*PZ+P.Z0*(1-PZ)


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2-1,beta=beta.old)
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

    #P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1-1,beta=beta.old))
    #PP<-P.Z1*PZ+P.Z0*(1-PZ)


    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1-1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff*Sweights.1
    #ff<-ff/PP*Sweights.1


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



EM.cc.CPV.Probit.Small.Nonpar2.EL.Weighted.X.Short.New<-function(Z,Sout,S1,W,Y,X,Wu,Xu,delta1,delta2,pd1.YZX,beta){

   #Sout=Sout.NC.All.X;pd1.YZX=pd1.YZX.true;pd2.YZXW=pd2.YZXW.true;beta=gamma;



      PZ=mean(Z)

        ## probability of delta_1=1 given Y=y,Z=z,X=x
       # pi.Y0Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,Xu)
#        pi.Y1Z0Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,Xu)
#        pi.Y0Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,Xu)
#        pi.Y1Z1Xu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,Xu)


        pi.Y0Z0Xu<-pd1.YZX[1,1,]
        pi.Y1Z0Xu<-pd1.YZX[2,1,]
        pi.Y0Z1Xu<-pd1.YZX[1,2,]
        pi.Y1Z1Xu<-pd1.YZX[2,2,]

        pi.XYZ<-c(pi.Y0Z0Xu,pi.Y1Z0Xu,pi.Y0Z1Xu,pi.Y1Z1Xu)


        WXu.ep<-expand.grid(Wu,Xu)


       #    pi.Y1Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,1)
#           pi.Y0Z1WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,1,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,1)
#           pi.Y1Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,1,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],1,0)
#           pi.Y0Z0WXu<-calPd1.ZYX.sub(Y,Z,X,delta1,0,0,WXu.ep[,2])*calPd2.d1YZWX.sub(W,Y,Z,X,delta1,delta2,WXu.ep[,1],WXu.ep[,2],0,0)

       pi.Y1Z1WXu<-pd2.YZXW[2,2,]
       pi.Y0Z1WXu<-pd2.YZXW[1,2,]
       pi.Y1Z0WXu<-pd2.YZXW[2,1,]
       pi.Y0Z0WXu<-pd2.YZXW[1,1,]



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



    oo<-match(WX[delta==0],WXu)

    ind<-rep(1:sum(delta==0), unlist(lval[oo]))


    YY<-rep(Y[delta==0],unlist(lval[oo]))


    WW<-rep(W[delta==0],unlist(lval[oo]))


    ZZ<-rep(Z[delta==0],unlist(lval[oo]))

    XX<-rep(X[delta==0],unlist(lval[oo]))


    SS1<-unlist(Sval[oo])

    Sweights<-unlist(Sout$Sweights[oo])

    WWXX<-WW+XX*10


   #### probability of P(Z=z|X=x)
   #PZ.X=table(X,Z)/as.numeric(table(X))

   oo.x=match(XX,Xu)

    beta.old<-beta



    oo.wx<-match(WWXX,WXu)
    pi.Y0Z0WWXX=pi.Y0Z0WXu[oo.wx]
    pi.Y1Z0WWXX=pi.Y1Z0WXu[oo.wx]
    pi.Y0Z1WWXX=pi.Y0Z1WXu[oo.wx]
    pi.Y1Z1WWXX=pi.Y1Z1WXu[oo.wx]



    XYZu<-expand.grid(Xu,c(0,1),c(0,1))
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
   # P.Z1<-pi.Y1Z1WWXX*R.probit(Z=1,S1=SS1,W=WW-1,beta=beta.old)+pi.Y0Z1WWXX*(1-R.probit(Z=1,S1=SS1,W=WW-1,beta=beta.old))


   # P.Z0<-pi.Y0Z0WWXX*(1-R.probit(Z=0,S1=SS1,W=WW-1,beta=beta.old))
   # PP<-P.Z1*PZ+P.Z0*(1-PZ)


    ff<-R.probit(Z=ZZ,S1=SS1,W=WW-1,beta=beta.old)
    ff<-ifelse(YY==1,ff,1-ff)
    ff<-ff*Sweights
    #ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind,sum)

     denom<-denom[ind]
     weights<-as.numeric(ff)/denom
     } else weights<-ff

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



EM.cc.CPV.Probit.Small.Nonpar2.ELMix.SubW.X.Short.New<-function(Z,Sout,Sout.EL,Wout.EL,S1,W,Y,X,Wu,Xu,delta1,delta2,pd2.YZXW,beta){

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
    dataLong<-getDataLong(Y[delta1==0],Z[delta1==0],X[delta1==0],Wu,Xu,Sout.EL,Wout.EL)


#> getDataLong
#function(Y,Z,X,Wu,Xu,Sout,Wout){
#
##Y=Y[delta1==0];Z=Z[delta1==0];X=X[delta1==0]
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


    #### probability of P(Z=z|X=x)
    #PZ.X=table(X,Z)/as.numeric(table(X))

    beta.old<-beta

    iter=0
    repeat{

    ### probability of delta2=1 conditional on S1,W,X

    P.Z1<-pd2.YZXW[2,2,oo.22]*R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old)+pd2.YZXW[1,2,oo.22]*(1-R.probit(Z=1,S1=SS1.2,W=WW.2-1,beta=beta.old))

    P.Z0<-pd2.YZXW[1,1,oo.22]*(1-R.probit(Z=0,S1=SS1.2,W=WW.2-1,beta=beta.old))


    PP<-P.Z1*PZ+P.Z0*(1-PZ)


    ff<-R.probit(Z=ZZ.2,S1=SS1.2,W=WW.2-1,beta=beta.old)
    ff<-ifelse(YY.2==1,ff,1-ff)
    ff<-ff/PP



     if (length(ff)>0){
     denom<-tapply(ff,ind.2,sum)
     denom<-denom[ind.2]
     weights.2<-as.numeric(ff)/denom
     } else weights.2<-ff


    ### probability of delta2=1 conditional on S1,W,X

    #P.Z1<-pd2.YZXW[2,2,oo.11]*R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old)+pd2.YZXW[1,2,oo.11]*(1-R.probit(Z=1,S1=SS1.1,W=WW.1-1,beta=beta.old))

    #P.Z0<-pd2.YZXW[1,1,oo.11]*(1-R.probit(Z=0,S1=SS1.1,W=WW.1-1,beta=beta.old))
    #PP<-P.Z1*PZ+P.Z0*(1-PZ)


    ff<-R.probit(Z=ZZ.1,S1=SS1.1,W=WW.1-1,beta=beta.old)
    ff<-ifelse(YY.1==1,ff,1-ff)
    ff<-ff*Sweights.1
    #ff<-ff/PP*PP.delta2*Sweights.1
    #ff<-ff/PP*Sweights.1


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
