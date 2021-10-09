### computer a few probabilties to enter into the pseudo-score estimation



############### Estimate p.delta2 (have to be estimated in practice since it won't be known) ######################3

   kk=table(delta2[delta1==1],X[delta1==1],W[delta1==1])


        p.delta2<-rep(NA,length(Xu)*length(Wu))
    for (j in 1:length(Xu)){
    for (k in 1:length(Wu)){
     p.delta2[(j-1)*length(Wu)+k]<-kk[2,j,k]/sum(kk[,j,k])
     #print((j-1)*length(Wu)+k)
    }
    }

  kk.long<-table(delta2[delta1==1],Y[delta1==1],Z[delta1==1],X[delta1==1],W[delta1==1])


 ############## Estimate pd2.YZXW ##############################

  ### pd2.d1YZXW is p(delta2=1|delta1=1,Y,Z,X,W)
    pd2.d1YZXW<-apply(kk.long,c(2,3,4,5),function(ci) {
    if (ci[1]+ci[2]>0) {
        ci[2]/(ci[1]+ci[2])} else {
      0
    }
    })



  # pd2.d1YZXW.true=pd2.d1YZXW
#   pd2.d1YZXW.true[1,1,,]<-ratio.p/ratio.p.w
#   pd2.d1YZXW.true[1,2,,]<-ratio.v/ratio.v.w


#   > pd2.d1YZXW[1,2,3,2]
#[1] 0.6142857
# sum(delta2==1 & delta1==1 & Y==0 & Z==1 & X==3 & W==2)/sum(delta1==1 & Y==0 & Z==1 & X==3 & W==2)

    kk.short<-table(delta1,Y,Z,X)
    kk.short.sum<-apply(kk.short,c(2,3,4),sum)


   # pd1.YZX is p(delta_1=1|Y,Z,X)

   if (length(unique(delta1))>1){
    pd1.YZX<-apply(kk.short,c(2,3,4),function(ci) ci[2]/(ci[1]+ci[2]))
    } else {
       pd1.YZX<-apply(kk.short,c(2,3,4),function(ci) 1)
    }


 #  if (length(unique(delta1))>1){
#      pd1.YZX<-apply(kk.short,c(2,3,4),function(ci) {
#
#         ifelse(ci[1]+ci[2]>0,ci[2]/(ci[1]+ci[2]),0)
#     }) } else {
#        pd1.YZX<-apply(kk.short,c(2,3,4),function(ci) 1)
#     }
#
#    pd1.YZX.true<-pd1.YZX

 #  for (k in 1:length(Xu))
#   pd1.YZX.true[,,k]<-matrix(c(ratio.p.w,ratio.v.w,1,1),byrow=T,nrow=2)

   for (i in 1:dim(pd1.YZX)[1]){
      for (j in 1:dim(pd1.YZX)[2]){
        for (k in 1:dim(pd1.YZX)[3]){
      if (is.na(pd1.YZX[i,j,k])) pd1.YZX[i,j,k]=0
      }
    }
    }  

   #> pd1.YZX[1,1,3]
   #[1] 0.4848485
   #> sum(delta1==1 & Y==0 & Z==0 & X==3)/sum(Y==0 & Z==0 & X==3)
   #[1] 0.4848485

   # pd2.YZXW= p(delta_2=1,delta_1=1|Y,Z,X,W)
    pd2.YZXW<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
    for (j in 1:length(Xu)){
    for (k in 1:length(Wu)){
     pd2.YZXW[,,(j-1)*length(Wu)+k]=pd2.d1YZXW[,,j,k]*pd1.YZX[,,j]
    }
    }


 #    pd2.YZXW.true<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
#    for (j in 1:length(Xu)){
#    for (k in 1:length(Wu)){
#     pd2.YZXW.true[,,(j-1)*length(Wu)+k]=pd2.d1YZXW.true[,,j,k]*pd1.YZX.true[,,j]
#    }
#    }


  #  pd2.YZXW.true1<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
#    for (j in 1:length(Xu)){
#    for (k in 1:length(Wu)){
#     pd2.YZXW.true1[,,(j-1)*length(Wu)+k]=pd2.d1YZXW[,,j,k]*pd1.YZX.true[,,j]
#    }
#    }
#
#       pd2.YZXW.true2<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
#    for (j in 1:length(Xu)){
#    for (k in 1:length(Wu)){
#     pd2.YZXW.true2[,,(j-1)*length(Wu)+k]=pd2.d1YZXW.true[,,j,k]*pd1.YZX[,,j]
#    }
#    }
#



 #pY.Z0XW.true<-c(
#  0.14335033, 0.12782308, 0.11835199, 0.10682200, 0.13768620, 0.12374996,
#  0.11494740, 0.10306712, 0.13318868, 0.12024566, 0.11128969, 0.09851838)
# pY.Z1XW.true=c(
#  0.08791712, 0.06510082, 0.05321513, 0.04050273, 0.07910427, 0.06020211,
#  0.04914115, 0.03684689, 0.07249290, 0.05569731, 0.04529818, 0.03273191)

  #     pZ.true=0.5
#    p.delta2.true<-rep(NA,length(Xu)*length(Wu))
#    for (j in 1:length(Xu)){
#    for (k in 1:length(Wu)){
#     oo=(j-1)*length(Wu)+k
#     num=pd2.YZXW.true[1,1,oo]*(1-pY.Z0XW.true[oo])*(1-pZ.true)+pd2.YZXW.true[2,1,oo]*pY.Z0XW.true[oo]*(1-pZ.true)+
#pd2.YZXW.true[1,2,oo]*(1-pY.Z1XW.true[oo])*(pZ.true)+pd2.YZXW.true[2,2,oo]*pY.Z1XW.true[oo]*(1-pZ.true)
#     denom=pd1.YZX.true[1,1,j]*(1-pY.Z0XW.true[oo])*(1-pZ.true)+pd1.YZX.true[2,1,j]*pY.Z0XW.true[oo]*(1-pZ.true)+
#pd1.YZX.true[1,2,j]*(1-pY.Z1XW.true[oo])*(pZ.true)+pd1.YZX.true[2,2,j]*pY.Z1XW.true[oo]*(1-pZ.true)
#
#     p.delta2.true[oo]<-num/denom
#     #print((j-1)*length(Wu)+k)
#    }
#    }


###
 PZ.X=table(X,Z)/as.numeric(table(X))

 #PZ.X.true=matrix(rep(c(1-pZ.true,pZ.true),length(Xu)),byrow=T,ncol=2)
