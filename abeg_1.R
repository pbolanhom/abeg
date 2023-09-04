
#Packages
#install.packages(deSolve)

library(deSolve)

#Model 
system = function(t,y,parms){
  with(as.list(c(y,parms)),
       {
         dC=C*(p*exp(-((x1-z)^2)/(2*s1^2))*a*R1+(1-p)*exp(-((x2-z)^2)/(2*s2^2))*a*R2-d)  
         dR1=R1*(r1-l*R1-C*a*p)
         dR2=R2*(r2-l*R2-C*a*(1-p))            
         dp=tau*p*(1-p)*a*(exp(-((x1-z)^2)/(2*s1^2))*R1-exp(-((x2-z)^2)/(2*s2^2))*R2)                       
         return(list(c(C=dC,R1=dR1,R2=dR2,p=dp)))
       })
}

parm=c(x1=-1,x2=1,s1=1,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,z=0.3) #parameters
attach(as.list(parm))
n0 = c(C=0.5,R1=1,R2=1,p=0.5) #initial conditions
out=ode(y=n0,times=seq(from=1,to=1500,by=0.1) ,func=system,parms=parm) #solve the model

#getting the numerical results after the dynamics reaches the equilibrium
Ct=tail(out[,2])[1]
R1t=tail(out[,3])[1]
R2t=tail(out[,4])[1]
pt=tail(out[,5])[1]

#The ecological stability of generalism (calculated points)
e1=exp(-((x1-z)^2)/(2*s1^2))
e2=exp(-((x2-z)^2)/(2*s2^2))
##The generalist scenario
R1G=d/(a*e1)
R2G=d/(a*e2)
pG=e2*(e1*a*r1-l*d)/(e1*e2*a*(r1+r2)-l*d*(e1+e2))
CG=(r1+r2)/a*(1-(l*d*(e1+e2))/(a*e1*e2*(r1+r2)))

##comparing the temporal and calculated points
iet=c("C"=Ct,"R1"=R1t,"R2"=R2t,"p"=pt)
ieG=c("C"=CG,"R1"=R1G,"R2"=R2G,"p"=pG)
data.frame(iet,ieG)

#Temporal dynamics
#Time Series
matplot(x=out[,1],y=out[,c(2:4)],type = "l",lwd=2,lty=c(1,2,2),col=c("black","darkgray"),
        xlab="Time",ylab="Population Density",frame=T)
legend("topright", legend = c("Consumer","Resource 1", "Resource 2"),lty=c(1,2,2), lwd=1, col=c("black","darkgray"),cex=0.7)
abline(h=c(CG,R1G,R2G),col="red",lty=2)
matplot(x=out[,1],y=out[,5],type = "l",lwd=3,lty=1,
        xlab="Time",ylab="Consumer foraging effort over Resource 1",frame=T,ylim=c(0,1))
abline(h=pG,col="red",lty=2)

#EVOLUTIONARY DYNAMICS
##The generalist singularity


#evaluating the zeta function
##defining the functions to be iterated
zeta0=function(z.initial) {
  e11=exp(-((parm[[1]]-z.initial)^2)/(2*parm[[3]]^2))
  e22=exp(-((parm[[2]]-z.initial)^2)/(2*parm[[4]]^2))
  p=e22*(e11*parm[[5]]*parm[[7]]-parm[[9]]*parm[[6]])/(e11*e22*parm[[5]]*(parm[[7]]+parm[[8]])-parm[[9]]*parm[[6]]*(e11+e22))
  (p*parm[[1]]*parm[[4]]^2+(1-p)*parm[[2]]*parm[[3]]^2)/(p*parm[[4]]^2+(1-p)*parm[[3]]^2) 
}
#the actual iteration of the zeta function
zeta=function(z.initial){
  for(i in 1:50) {
    prev <- z.initial
    z.initial <- zeta0(z.initial)
    cat(i, z.initial, "\n")
    if (abs(z.initial - prev) < 1e-5) break
  }
  z0=z.initial
  return(z0)}
z.initial=parm[[11]] #setting the initial z value
z0=zeta(z.initial)

#######################################################################
#Proving that convergence leads to the same trait value as the eco-evolutionary model
system2 = function(t,y,parms){
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
parm2=c(x1=-1,x2=1,s1=1,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2,V=0.05) #have the same parameters values as the previous model
n02 = c(C=0.5,R1=1,R2=1,p=0.5,z=0.3)
time=seq(from=1,to=1500,by=0.1)
out2=ode(y=n02,times=time,func=system2,parms=parm2)
zt=tail(out2[,6])[1]
data.frame(zt,z0) #when the attractor is zero the values may not be identical. However, both are zero at the limit 



