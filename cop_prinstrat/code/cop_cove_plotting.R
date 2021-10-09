##### 

lwd=2;lwd2=2

for (a in assays){

ps<-c("Day"%.%pop%.%a)

load(file=paste0(save.results.to,"outCOVE_boot_",ps,"_",y.v,".Rdata"))

mypdf(file=paste0(save.results.to,"fig_",ps,'_',y.v,'_',study_name),width=6,height=6)
#plot(outb$Su,outb$VE.S,type='l',ylim=c(-1,1),xlab=quote(s[1]),ylab=quote(VE(s[1])),main=paste0(ps,',',y.v),xaxt="n", yaxt="n", col="red3", lwd=2.5)
plot(outb$Su,outb$VE.S,type='l',ylim=c(0.5,1),xlab="S(1)",ylab="VE(S(1))",main=paste0(ps,',',y.v),xaxt="n", yaxt="n", col="red3", lwd=2.5)
lines(outb$Su,outb$VE.S,col=2,lwd=lwd)
lines(outb$Su,outb$low.ci,lty=2,lwd=lwd)
lines(outb$Su,outb$high.ci,lty=2,lwd=lwd)

lines(outb$Su,outb$low.cb,lty=3,lwd=lwd)
lines(outb$Su,outb$high.cb,lty=3,lwd=lwd2)
axis(side=1, at=seq(0.5,5.5,by=0.5))
axis(side=2, at=seq(-1,1,by=0.25), labels=paste0(seq(-100,100,by=25),"%"))
#mtext("Estimated VE by Titer", side=2, las=0, line=3, cex=1.1)
#mtext("VE(s1)", side=2, las=0, line=3, cex=1.1)

#legend('topleft', fill=colHist, border=colHist, legend="Vaccine Group", bty="n", cex=1)
legend(min(outb$Su),-0.5,lty=2,cex=1,c('95% pointwise CI'),bty='n',lwd=lwd)
legend(min(outb$Su),-0.7,lty=3,cex=1,c('95% simultaneous CI'),bty='n',lwd=lwd2)
legend('bottom',paste0('p.int=',round(outb$pvalue,3)),bty='n')
#legend(2.8,1.9,paste('Hinge point=',round(thre,2),sep=''),bty='n',cex=1.1)
dev.off()

}
