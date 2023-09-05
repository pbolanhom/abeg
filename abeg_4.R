
#####################################################################################################
###pairwise-invasibility analysis
#Model #with trait evolution (z is a parameter)
system0 = function(t,y,parms){
  with(as.list(c(y,parms)),
       {
         dC=C*(p*exp(-((x1-z)^2)/(2*s1^2))*a*R1+(1-p)*exp(-((x2-z)^2)/(2*s2^2))*a*R2-d)  
         dR1=R1*(r1-l*R1-C*a*p)
         dR2=R2*(r2-l*R2-C*a*(1-p))            
         dp=tau*p*(1-p)*a*(exp(-((x1-z)^2)/(2*s1^2))*R1-exp(-((x2-z)^2)/(2*s2^2))*R2)                      
         return(list(c(C=dC,R1=dR1,R2=dR2,p=dp)))
       })
}
parm0=c(x1=-1,x2=1,s1=0.93,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2)
detach(as.list(parm0))
attach(as.list(parm0))

zz = seq(from = -2, to=2, by=0.1)
zCminmax = matrix(NA, ncol=2, nrow=length(zz))
zR1minmax = matrix(NA, ncol=2, nrow=length(zz))
zR2minmax = matrix(NA, ncol=2, nrow=length(zz))
zpminmax = matrix(NA, ncol=2, nrow=length(zz))
for(i in 1:length(zz)){
  parmsi=c(x1=-1,x2=1,s1=1,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,z=zz[i])
  n = c(C=0.5,R1=1,R2=1,p=0.5)
  outi=ode(y=n,times=seq(from=1,to=1000,by=0.1),func=system0,parms=parmsi)
  outi.log = outi#log(outi+1)
  zCminmax[i,] = range(outi.log[(nrow(outi.log)-500):nrow(outi.log),2])
  zR1minmax[i,] = range(outi.log[(nrow(outi.log)-500):nrow(outi.log),3])
  zR2minmax[i,] = range(outi.log[(nrow(outi.log)-500):nrow(outi.log),4])
  zpminmax[i,] = range(outi.log[(nrow(outi.log)-500):nrow(outi.log),5])
}

e1MUT=exp(-((x1-zz)^2)/(2*s1^2))
e2MUT=exp(-((x2-zz)^2)/(2*s2^2))
pres=matrix(zpminmax[,1])
R1res=matrix(zR1minmax[,1])
R2res=matrix(zR2minmax[,1])
Fi=NULL
Fi0=NULL
for (i in 1:length(pres)) {
  Fi0=pres[i]*e1MUT*a*R1res[i]+(1-pres[i])*e2MUT*a*R2res[i]-d
  Fi=cbind(Fi,Fi0)
}
#FOR THE ESS CASE
invasPIP11=as.data.frame(Fi)
invasPIP=as.matrix(invasPIP0)
colnames(invasPIP11)=as.numeric(zz)
rownames(invasPIP11)=as.numeric(zz)
PIP11=invasPIP11>0
pip11=PIP11*1
#library(plot.matrix)
#library('psych')
head(pip11)
par()
#plot(pip,border=NA,asp=T,col=c("white","black"),main="",key=NULL)
plot(1, type = "n", xlab = "",ylab="",xlim=c(-2, 2), ylim = c(-2, 2))
image(zz,zz,t(pip11),col=c("white","black"),useRaster = 1,add=T)
abline(v=0,lty=2,lwd=2,col="darkgray")


#FOR THE BP CASE
#invasPIP22=as.data.frame(Fi)
#invasPIP=as.matrix(invasPIP0)
colnames(invasPIP22)=as.numeric(zz)
rownames(invasPIP22)=as.numeric(zz)
PIP22=invasPIP22>0
pip22=PIP22*1
#library(plot.matrix)
#library('psych')
head(pip22)
par()
#plot(pip,border=NA,asp=T,col=c("white","black"),main="",key=NULL)
plot(1, type = "n", xlab = "",ylab="",xlim=c(-2, 2), ylim = c(-2, 2))
image(zz,zz,t(pip22),col=c("white","black"),useRaster = 1,add=T)
abline(v=0,lty=2,lwd=2,col="darkgray")

#FOR THE repellor CASE
#invasPIP33=as.data.frame(Fi)
#invasPIP=as.matrix(invasPIP0)
colnames(invasPIP33)=as.numeric(zz)
rownames(invasPIP33)=as.numeric(zz)
PIP33=invasPIP33>0
pip33=PIP33*1
#library(plot.matrix)
#library('psych')
head(pip33)
par()
#plot(pip,border=NA,asp=T,col=c("white","black"),main="",key=NULL)
plot(1, type = "n", xlab = "",ylab="",xlim=c(-2, 2), ylim = c(-2, 2))
image(zz,zz,t(pip33),col=c("white","black"),useRaster = 1,add=T)
abline(v=0,lty=2,lwd=2,col="darkgray")

############
#fitness landscapes 

parm=c(x1=-1,x2=1,s1=1.1,s2=1.1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=0)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
ESSSIM=SIMIE(out)
attach(ESSSIM)
ESSSIM_ep=evopoints(R1t,R2t,pt,zt)
ESSSIM_ec=evocurves(R1t,R2t,pt)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=ESSSIM_ec$FitSurf,type = "l",lwd=3,lty=1,xlab="z",ylab="Fitness",
        frame=F,ylim=c(-0.1,0.1),xlim=c(-2,2))
matpoints(x=ESSSIM$zt,y=ESSSIM_ep$Fit,pch=19,col="black")
#abline(v=ESSSIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

#for bp
parm=c(x1=-1,x2=1,s1=0.93,s2=0.93,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=0)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
BPSIM=SIMIE(out)
attach(BPSIM)
BPSIM_ep=evopoints(R1t,R2t,pt,zt)
BPSIM_ec=evocurves(R1t,R2t,pt)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=BPSIM_ec$FitSurf,type = "l",lwd=3,lty=1,xlab="z",ylab="Fitness",
        frame=F,ylim=c(-0.1,0.1),xlim=c(-2,2))
matpoints(x=BPSIM$zt,y=BPSIM_ep$Fit,pch=19,col="black")
#abline(v=BPSIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

#for repellor
parm=c(x1=-1,x2=1,s1=0.91,s2=0.91,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=0)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
RPSIM=SIMIE(out)
attach(RPSIM)
RPSIM_ep=evopoints(R1t,R2t,pt,zt)
RPSIM_ec=evocurves(R1t,R2t,pt)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=RPSIM_ec$FitSurf,type = "l",lwd=3,lty=1,xlab="z",ylab="Fitness",
        frame=F,ylim=c(-0.1,0.1),xlim=c(-2,2))
matpoints(x=RPSIM$zt,y=RPSIM_ep$Fit,pch=19,col="black")
#abline(v=RPSIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

