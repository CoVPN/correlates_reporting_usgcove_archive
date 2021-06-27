### computer a few probabilties to enter into the pseudo-score estimation


    kk.short<-table(delta2,Y,Z,X)
    kk.short.sum<-apply(kk.short,c(2,3,4),sum)


   # pd2.YZX is p(delta_2=1|Y,Z,X)

   if (length(unique(delta2))>1){
    pd2.YZX<-apply(kk.short,c(2,3,4),function(ci) ci[2]/(ci[1]+ci[2]))
    } else {
       pd2.YZX<-apply(kk.short,c(2,3,4),function(ci) 1)
    }
 
 
   ## there will be NA if there does not exist certainly Y,Z,X combination in the data, just change it to zero
   for (i in 1:dim(pd2.YZX)[1]){
      for (j in 1:dim(pd2.YZX)[2]){
        for (k in 1:dim(pd2.YZX)[3]){
      if (is.na(pd2.YZX[i,j,k])) pd2.YZX[i,j,k]=0
      }
    }
    }  


 
    pd2.YZXW<-array(NA,dim=c(2,2,length(Xu)*length(Wu)))
    for (j in 1:length(Xu)){
    for (k in 1:length(Wu)){
     pd2.YZXW[,,(j-1)*length(Wu)+k]=pd2.YZX[,,j]
    }
    }
