rm(list=ls())
library(kyotil)

#index1=1;index2=1

setwd("C:/All_Files/1yingsstuff/correlates_reporting/cop_prinstrat")

assay.seq<-c("bindSpike","bindRBD","bindN","pseudoneutid50","pseudoneutid80","liveneutmn50")
pop.seq<-c('57','29')

ps.seq<-c('EIA.log10d14overd0','PCA.log10d14overd0','RSVA.log10d14','RSVB.log10d14')
bsm.seq<-c('EIA.log10d0','PCA.log10d0','RSVA.log10d0','RSVB.log10d0')
y.seq<-c('y1','y2')


ps.name<-c("bindSpike","bindRBD","bindN","pseudoneutid50","pseudoneutid80","liveneutmn50")
y.name<-c("Day 57", "Day 29")


getcicb<-function(Su,VE.out,VE.S,beta.out,beta){

Q.95.beta<-quantile(beta.out,0.95,na.rm=T)
Q.05.beta<-quantile(beta.out,.05,na.rm=T)

Q.9.beta<-quantile(beta.out,0.9,na.rm=T)
Q.1.beta<-quantile(beta.out,.1,na.rm=T)

#SE.beta<-(Q.95.beta-Q.05.beta)/qnorm(0.95)/2
SE.beta<-(Q.95.beta-Q.05.beta)/qnorm(0.95)/2

pvalue<-2*(1-pnorm(abs(beta/SE.beta)))

SE.RR=apply(log(1-VE.out),2,sd,na.rm=T)

Q.95=apply(log(1-VE.out),2,quantile,0.95,na.rm=T)
Q.05=apply(log(1-VE.out),2,quantile,0.05,na.rm=T)

Q.9=apply(log(1-VE.out),2,quantile,0.9,na.rm=T)
Q.1=apply(log(1-VE.out),2,quantile,0.1,na.rm=T)

SE.RR1=(Q.95-Q.05)/2/qnorm(0.95)
SE.RR2=(Q.9-Q.1)/2/qnorm(0.9)

U=rep(NA,nrow(VE.out))
for (kk in 1:nrow(VE.out)){
U[kk]=max(abs(log(1-VE.out[kk,])-log(1-VE.S))/SE.RR1)
}
c.alpha=quantile(U,0.95,na.rm=T)


high.ci=1-exp(log(1-VE.S)-qnorm(0.975)*SE.RR1)
low.ci=1-exp(log(1-VE.S)+qnorm(0.975)*SE.RR1)

high.cb=1-exp(log(1-VE.S)-c.alpha*SE.RR1)
low.cb=1-exp(log(1-VE.S)+c.alpha*SE.RR1)


return(list(Su=Su,VE.S=VE.S,low.ci=low.ci,high.ci=high.ci,low.cb=low.cb,high.cb=high.cb,pvalue=pvalue,SE.beta=SE.beta,SE.RR=SE.RR1))
}


lwd=2;lwd2=2


##### 

for (index1 in 1:2){
for (index2 in 1:2){


pop=pop.seq[index1]
y.v<-"EventIndPrimaryD"%.%pop       

ps<-c("Day"%.%pop%.%assay.seq[index2])

ps<-ps.seq[index1]
y.v<-y.seq[index2]

load(file=paste0("input/Result/outCohort_",ps,"_",y.v,".Rdata"))
load(file=paste0("input/Result/outbootCohort_",ps,"_",y.v,".Rdata"))

out=getcicb(VE1$Su,out.VE1.seq,VE1$VE,out.beta1.seq[,4],fit1[4])
low.ci=out$low.ci;high.ci=out$high.ci;Su=out$Su;VE.S=out$VE.S;low.cb=out$low.cb;high.cb=out$high.cb
mypdf(file=paste0("input/figCohort_nocov",ps.seq[index1],'_',y.seq[index2]),width=6,height=6)
plot(Su,VE.S,type='l',ylim=c(-1,1),xlab=quote(s[1]),ylab=quote(VE(s[1])),main=paste0(ps.name[index1],',',y.name[index2]),xaxt="n", yaxt="n", col="red3", lwd=2.5)
lines(Su,VE.S,col=2,lwd=lwd)
lines(Su,low.ci,lty=2,lwd=lwd)
lines(Su,high.ci,lty=2,lwd=lwd)

lines(Su,low.cb,lty=3,lwd=lwd)
lines(Su,high.cb,lty=3,lwd=lwd2)
axis(side=1, at=seq(0.5,5.5,by=0.5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
#mtext("Estimated VE by Titer", side=2, las=0, line=3, cex=1.1)
#mtext("VE(s1)", side=2, las=0, line=3, cex=1.1)

#legend('topleft', fill=colHist, border=colHist, legend="Vaccine Group", bty="n", cex=1)
legend(min(Su),-0.5,lty=2,cex=1,c('95% pointwise CI'),bty='n',lwd=lwd)
legend(min(Su),-0.7,lty=3,cex=1,c('95% simultaneous CI'),bty='n',lwd=lwd2)
legend('bottom',paste0('W/o cov adj, p.int=',round(out$pvalue,3)),bty='n')
#legend(2.8,1.9,paste('Hinge point=',round(thre,2),sep=''),bty='n',cex=1.1)
dev.off()

rm(low.ci,high.ci,low.cb,high.cb,VE.S,out)
 
out=getcicb(VE3$Su,out.VE3.seq,VE3$VE,out.beta3.seq[,4],fit3[4])
low.ci=out$low.ci;high.ci=out$high.ci;Su=out$Su;VE.S=out$VE.S;low.cb=out$low.cb;high.cb=out$high.cb

mypdf(file=paste0("input/figCohort_withcov",ps.seq[index1],'_',y.seq[index2]),width=6,height=6)
plot(Su,VE.S,type='l',ylim=c(-1,1),xlab=quote(s[1]),ylab=quote(VE(s[1])),main=paste0(ps.name[index1],',',y.name[index2]),xaxt="n", yaxt="n", col="red3", lwd=2.5)
lines(Su,VE.S,col=2,lwd=lwd)
lines(Su,low.ci,lty=2,lwd=lwd)
lines(Su,high.ci,lty=2,lwd=lwd)

lines(Su,low.cb,lty=3,lwd=lwd)
lines(Su,high.cb,lty=3,lwd=lwd2)
axis(side=1, at=seq(0.5,5.5,by=0.5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
#mtext("Estimated VE by Titer", side=2, las=0, line=3, cex=1.1)
#mtext("VE(s1)", side=2, las=0, line=3, cex=1.1)

#legend('topleft', fill=colHist, border=colHist, legend="Vaccine Group", bty="n", cex=1)
legend(min(Su),-0.5,lty=2,cex=1,c('95% pointwise CI'),bty='n',lwd=lwd)
legend(min(Su),-0.7,lty=3,cex=1,c('95% simultaneous CI'),bty='n',lwd=lwd2)
legend('bottom',paste0('With cov adj, p.int=',round(out$pvalue,3)),bty='n')
#legend(2.8,1.9,paste('Hinge point=',round(thre,2),sep=''),bty='n',cex=1.1)
dev.off()

rm(Su,out.beta1.seq,out.beta2.seq,out.beta3.seq,out.beta4.seq,
out.VE1.seq,out.VE2.seq,out.VE3.seq,out.VE4.seq,low.ci,high.ci,low.cb,high.cb,VE.S,out)
 
}
}
