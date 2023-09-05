#The adaptive dynamics of a behaviorally flexible generalist
#Packages
library(deSolve)
#Model #with trait evolution (z is a variable)
system = function(t,y,parms){
  with(as.list(c(y,parms)),
       {
         dC=C*(p*exp(-((x1-z)^2)/(2*s1^2))*a*R1+(1-p)*exp(-((x2-z)^2)/(2*s2^2))*a*R2-d)  
         dR1=R1*(r1-l*R1-C*a*p)
         dR2=R2*(r2-l*R2-C*a*(1-p))            
         dp=tau*p*(1-p)*a*(exp(-((x1-z)^2)/(2*s1^2))*R1-exp(-((x2-z)^2)/(2*s2^2))*R2)                      
         dz=V*a*(p*(x1-z)/(s1^2)*exp(-((x1-z)^2)/(2*s1^2))*R1+(1-p)*(x2-z)/(s2^2)*exp(-((x2-z)^2)/(2*s2^2))*R2)          
         return(list(c(C=dC,R1=dR1,R2=dR2,p=dp,z=dz)))
       })
}
parm=c(x1=-1,x2=1,s1=1,s2=1,a=1,d=0.4,r1=0.8,r2=1,l=0.25,tau=0,V=0.05)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=1)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
###########Functions
###Simulated interior equilibrium
SIMIE=function(OUT){
  Ct=tail(OUT[,2])[1]
  R1t=tail(OUT[,3])[1]
  R2t=tail(OUT[,4])[1]
  pt=tail(OUT[,5])[1]
  zt=tail(OUT[,6])[1]
  return(data.frame(Ct,R1t,R2t,pt,zt))
}
###Calculated interior equilibrium
CALCIE=function(Z){
  e1=exp(-((x1-Z)^2)/(2*s1^2))
  e2=exp(-((x2-Z)^2)/(2*s2^2))
  R1c=d/(a*e1)
  R2c=d/(a*e2)
  pc=e2*(e1*a*r1-l*d)/(e1*e2*a*(r1+r2)-l*d*(e1+e2))
  Cc=(r1*a*e1-l*d)/(e1*pt*a^2)
  ct1=(R1t*a*e1)/(s1^2)
  ct2=(R2t*a*e2)/(s2^2)
  zc=(pt*x1*ct1+(1-pt)*x2*ct2)/(pt*ct1+(1-pt)*ct2)
  return(data.frame(Cc,R1c,R2c,pc,zc))
}
###Finding generalist singularity
zStar=function(Z) {
  e11=exp(-((x1-Z)^2)/(2*s1^2))
  e22=exp(-((x2-Z)^2)/(2*s2^2))
  p=e22*(e11*a*r1-l*d)/(e11*e22*a*(r1+r2)-l*d*(e11+e22))
  (p*x1*s2^2+(1-p)*x2*s1^2)/(p*s2^2+(1-p)*s1^2) 
}
zeta=function(Z){
  for(i in 1:50) {
    z0=Z
    Z=zStar(Z)
    cat(i, Z, "\n")
    if (abs(Z - z0) < 1e-5) break
  }
  zgen=Z
  e1gen=exp(-((x1-zgen)^2)/(2*s1^2))
  e2gen=exp(-((x2-zgen)^2)/(2*s2^2))
  pgen=e2gen*(e1gen*a*r1-l*d)/(e1gen*e2gen*a*(r1+r2)-l*d*(e1gen+e2gen))
  R1gen=d/(a*e1gen)
  R2gen=d/(a*e2gen)  
  return(data.frame(zgen,pgen,R1gen,R2gen))
}  
###Values of evolutionary derivatives (EVOPOINTS)
evopoints=function(R1,R2,p,z){
  e10=exp(-((x1-z)^2)/(2*s1^2))
  e20=exp(-((x2-z)^2)/(2*s2^2))
  Fit=p*e10*a*R1+(1-p)*e20*a*R2-d
  gSel=p*(x1-z)/(s1^2)*a*e10*R1+(1-p)*(x2-z)/(s2^2)*a*e20*R2
  ESS=a*(e10*p*R1*(((x1-z)^2)/(s1^4)-1/(s1^2))+e20*(1-p)*R2*(((x2-z)^2)/(s2^4)-1/(s2^2)))
  num=(z-x2)*e20/(s2^2)+(z-x1)*e10/(s1^2)
  dR1p.dz=(d*e10*(z-x1))/(a*(d*l*(e20+e10)-r2-r1))*((d*l+(d*l*e10-a*r1))/(s1^2)-(d*l*(d*l*e10-a*r1))*(num)/((z-x1)*(d*l*(e20+e10)-r2-r1)))
  dR2p.dz=(d*e20*(z-x2))/(a)*(((d*l+(d*l*e20-r2+(a-1)*r1))/((d*l*(e20+e10)-r2-r1)*s2^2)-(d*l*(d*l*e20-r2+(a-1)*r1)))*(num)/((z-x2)*(d*l*(e20+e10)-r2-r1)^2))
  CSS0=a*((x1-z)/(s1^2)*e10*dR1p.dz+(x2-z)/(s2^2)*e20*dR2p.dz)
  CSS=CSS0+ESS
  return(data.frame(Fit,gSel,ESS,CSS))
} 
###Invasion fitness surfaces (EVOCURVES)
evocurves=function(R1,R2,p){
  zz = seq(from = -2, to=2, by=0.01)
  e1zz=exp(-((x1-zz)^2)/(2*s1^2))
  e2zz=exp(-((x2-zz)^2)/(2*s2^2))
  FitSurf=p*e1zz*a*R1+(1-p)*e2zz*a*R2-d
  FitSurfSP1=e1zz*a*R1-d
  FitSurfSP2=e2zz*a*R2-d
  SGzz=p*(x1-zz)/(s1^2)*a*e1zz*R1+(1-p)*(x2-zz)/(s2^2)*a*e2zz*R2
  return(data.frame(FitSurf,FitSurfSP1,FitSurfSP2,SGzz))
}



#Fig 1
#with generalist convergence
#V=0
parm=c(x1=-1,x2=1,s1=0.8,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=-1)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
cV0SIM=SIMIE(out)
attach(cV0SIM)
cV0SIM_ep=evopoints(R1t,R2t,pt,zt)
cV0SIM_ec=evocurves(R1t,R2t,pt)

#V=0.05
parm=c(x1=-1,x2=1,s1=0.8,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0.05)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=-1)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
cV05SIM=SIMIE(out)
attach(cV05SIM)
cV05SIM_ep=evopoints(R1t,R2t,pt,zt)
cV05SIM_ec=evocurves(R1t,R2t,pt)


#NO CONVERGENCE
#V=0
parm=c(x1=-1.5,x2=1,s1=0.8,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=-1.5)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
dV0SIM=SIMIE(out)
attach(dV0SIM)
dV0SIM_ep=evopoints(R1t,R2t,pt,zt)
dV0SIM_ec=evocurves(R1t,R2t,pt)
#V=0.05
parm=c(x1=-1.5,x2=1,s1=0.8,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0.05)
detach(as.list(parm))
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5,z=-1.5)
time=seq(from=1,to=1500,by=0.1)
out=ode(y=n0,times=time,func=system,parms=parm)
dV05SIM=SIMIE(out)
attach(dV05SIM)
dV05SIM_ep=evopoints(R1t,R2t,pt,zt)
dV05SIM_ec=evocurves(R1t,R2t,pt)




zz = seq(from = -2, to=2, by=0.01)
#Fig fitness a)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=cV05SIM_ec$FitSurf,type = "l",lwd=3,lty=1,xlab="z",ylab="Fitness",
        frame=T,ylim=c(-0.1,0.4),xlim=c(-2,2))
lines(x=zz,y=cV0SIM_ec$FitSurf,type = "l",lwd=3,lty=1,col="red")
matpoints(x=cV05SIM$zt,y=cV05SIM_ep$Fit,pch=19,col="black")
abline(v=cV0SIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

#fig grad a)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=cV05SIM_ec$SGzz,type = "l",lwd=3,lty=1,xlab="z",ylab="Selection gradient",
        frame=F,ylim=c(-0.6,0.6),xlim=c(-2,2))
lines(x=zz,y=cV0SIM_ec$SGzz,type = "l",lwd=3,lty=1,col="darkgray")
matpoints(x=cV05SIM$zt,y=cV05SIM_ep$gSel,pch=19,col="black")
abline(v=cV0SIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

#Fig fitness b)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=dV05SIM_ec$FitSurf,type = "l",lwd=3,lty=1,xlab="z",ylab="Fitness",
        frame=T,ylim=c(-0.1,0.4),xlim=c(-2,2))
lines(x=zz,y=dV0SIM_ec$FitSurf,type = "l",lwd=3,lty=1,col="red")
matpoints(x=dV05SIM$zt,y=dV05SIM_ep$Fit,pch=17,col="black")
#matpoints(x=dV0SIM$zt,y=dV0SIM_ep$Fit,pch=19,col="black")
abline(v=dV0SIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")
#fig grad b)
par(mar=c(5,5,1,1)+0.1,mgp=c(3,1,0),cex=2,xpd=FALSE,las=1)
matplot(x=zz,y=dV05SIM_ec$SGzz,type = "l",lwd=3,lty=1,xlab="z",ylab="Selection gradient",
        frame=F,ylim=c(-0.6,0.6),xlim=c(-2,2))
lines(x=zz,y=dV0SIM_ec$SGzz,type = "l",lwd=3,lty=1,col="darkgray")
matpoints(x=dV05SIM$zt,y=dV05SIM_ep$gSel,pch=19,col="black")
abline(v=dV0SIM$zt,lty=2,col="black",lwd=2)
abline(h=0,lty=1,col="black")

