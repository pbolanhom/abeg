

#The evo-behavior feedback loop
#Here, we show that behavior update function zeta produces an evo-behavioral feedback loop that 
# generates the same interior equilibrium as an eco-evolutionary system 
#(this is only valie for cases in which the generalist singularity exists). 
#This means that, in our model,
#the evo-behavioral feedback dampens the eco-evolutionary loop. 


#Packages
#install.packages(deSolve)
library(deSolve)
#No eco-evo loop
system.1 = function(t,y,parms){
  with(as.list(c(y,parms)),
       {
         dC=C*(p*exp(-((x1-z)^2)/(2*s1^2))*a*R1+(1-p)*exp(-((x2-z)^2)/(2*s2^2))*a*R2-d)  
         dR1=R1*(r1-l*R1-C*a*p)
         dR2=R2*(r2-l*R2-C*a*(1-p))            
         dp=tau*p*(1-p)*a*(exp(-((x1-z)^2)/(2*s1^2))*R1-exp(-((x2-z)^2)/(2*s2^2))*R2)                       
         return(list(c(C=dC,R1=dR1,R2=dR2,p=dp)))
       })
}

#With eco-evo loop
system.2 = function(t,y,parms){
  with(as.list(c(y,parms)),
       {
         dC=C*(p*exp(-((x1-Z)^2)/(2*s1^2))*a*R1+(1-p)*exp(-((x2-Z)^2)/(2*s2^2))*a*R2-d)  
         dR1=R1*(r1-l*R1-C*a*p)
         dR2=R2*(r2-l*R2-C*a*(1-p))            
         dp=tau*p*(1-p)*a*(exp(-((x1-Z)^2)/(2*s1^2))*R1-exp(-((x2-Z)^2)/(2*s2^2))*R2)                      
         dZ=V*a*(p*(x1-Z)/(s1^2)*exp(-((x1-Z)^2)/(2*s1^2))*R1+(1-p)*(x2-Z)/(s2^2)*exp(-((x2-Z)^2)/(2*s2^2))*R2)          
         return(list(c(C=dC,R1=dR1,R2=dR2,p=dp,Z=dZ)))
       })
}

parameters=c(x1=-1,x2=1,s1=1,s2=1,a=1,d=0.4,r1=1,r2=1,l=0.25,tau=2)
attach(as.list(parameters))
V=0.05 #for the eco-evo system (system.2)
z=0.5 #for the eco system (system.1)
n01 = c(C=0.5,R1=1,R2=1,p=0.5)
n02 = c(C=0.5,R1=1,R2=1,p=0.5,Z=0.3)
time=seq(from=1,to=1500,by=0.1)
out1=ode(y=n01,times=time,func=system.1,parms=parameters)
out2=ode(y=n02,times=time,func=system.2,parms=parameters)

#taking the equilibrium solutions
#for sys 1
Ct1=tail(out1[,2])[1]
R1t1=tail(out1[,3])[1]
R2t1=tail(out1[,4])[1]
pt1=tail(out1[,5])[1]
iet1=c("C"=Ct1,"R1"=R1t1,"R2"=R2t1,"p"=pt1,"z"=z)
#for sys 2
Ct2=tail(out2[,2])[1]
R1t2=tail(out2[,3])[1]
R2t2=tail(out2[,4])[1]
pt2=tail(out2[,5])[1]
Zt2=tail(out2[,6])[1]
iet2=c("C"=Ct2,"R1"=R1t2,"R2"=R2t2,"p"=pt2,"z"=Zt2)
#evaluating the zeta function (for system 1)
##defining the functions to be iterated
zeta0=function(z.initial) {
  e11=exp(-((parameters[[1]]-z.initial)^2)/(2*parameters[[3]]^2))
  e22=exp(-((parameters[[2]]-z.initial)^2)/(2*parameters[[4]]^2))
  p=e22*(e11*parameters[[5]]*parameters[[7]]-parameters[[9]]*parameters[[6]])/(e11*e22*parameters[[5]]*(parameters[[7]]+parameters[[8]])-parameters[[9]]*parameters[[6]]*(e11+e22))
  (p*parameters[[1]]*parameters[[4]]^2+(1-p)*parameters[[2]]*parameters[[3]]^2)/(p*parameters[[4]]^2+(1-p)*parameters[[3]]^2) 
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
z0=zeta(z)

#The ecological stability of generalism (calculated points)
e1z0=exp(-((x1-z0)^2)/(2*s1^2))
e2z0=exp(-((x2-z0)^2)/(2*s2^2))
##The generalist scenario
R1G0=d/(a*e1z0)
R2G0=d/(a*e2z0)
pG0=e2z0*(e1z0*a*r1-l*d)/(e1z0*e2z0*a*(r1+r2)-l*d*(e1z0+e2z0))
CG0=(r1+r2)/a*(1-(l*d*(e1z0+e2z0))/(a*e1z0*e2z0*(r1+r2)))
iet1G0=c("C"=CG0,"R1"=R1G0,"R2"=R2G0,"p"=pG0,"z"=z0)
db=data.frame("sys1"=iet1,"sys2"=iet2,"sys1update"=iet1G0)
db
#Note that, in the db object, the second and third columns have the same values (when the generalist singularity exists). 
#small differences are due to rounding discrepancies between the numerical simulations and calculations.


