#test file (reduced model: Ken paper)
options(max.print=500)
HUB=read.csv("HUB.csv",header=T)

table(HUB$seq)

#R codes for generating results in Table 1
#estimate Q2-function (saturate model)  **checked**
fit2=lm(Y~A1+R+A2+A1*R*A2,data=HUB) 
B=coef(fit2)                                                           #theta2_hat
round(B,2)
SE2=summary(fit2)$sigma                                                #sigma2_hat

#estimate the value of Q2 function given (saturated model)  **checked**
Q2=function(a1,r,a2,beta){
   beta[which(is.na(beta))]=0 
   EY=beta[1]+beta[2]*a1+beta[3]*r+beta[4]*a2+
      beta[5]*a1*r+beta[6]*a1*a2+B[7]*r*a2+
      beta[8]*a1*r*a2
   return(as.numeric(EY))
}
Q2(a1=0,r=0,a2=0,beta=B)
Q2(a1=0,r=0,a2=1,beta=B)
Q2(a1=0,r=1,a2=0,beta=B)
Q2(a1=0,r=1,a2=1,beta=B)
Q2(a1=1,r=0,a2=0,beta=B)
Q2(a1=1,r=0,a2=1,beta=B)
Q2(a1=1,r=1,a2=0,beta=B)
Q2(a1=1,r=1,a2=1,beta=B)

#function to calculate pseudo maximum value of Q2 (pick the better one)  **checked**
Q2max=function(GA1,GR,beta){
      Q20=Q2(a1=GA1,r=GR,a2=0,b=beta)            #Q2 for A2=0 given a1 and r
      Q21=Q2(a1=GA1,r=GR,a2=1,b=beta)            #Q2 for A2=1 given a1 and r
      return(max(Q20,Q21))                       #returne the max value of Q2
}
Q2max(GA1=0,GR=0,beta=B)
Q2max(GA1=0,GR=1,beta=B)
Q2max(GA1=1,GR=0,beta=B)
Q2max(GA1=1,GR=1,beta=B)


#calcuate pseudo E(Y|A1,R) for estimating Q1-function
for (k in 1:nrow(HUB)){HUB$EQ2max[k]=Q2max(GA1=HUB$A1[k],GR=HUB$R[k],beta=B)}

#estimate Q1-function
fit1=lm(EQ2max~A1,data=HUB)                        #estiamte stage-1 Q-function
G=coef(fit1)                                       #theta1_hat
SE1=summary(fit1)$sigma                            #sigma t=1

#estimte value of Q1-function given A1             **checked**
Q1=function(a1,gamma){
   gamma[which(is.na(gamma))]=0 
   EY=gamma[1]+gamma[2]*a1
   return(as.numeric(EY))
}
Q1(a1=0,gamma=G)
Q1(a1=1,gamma=G)


#estimate optimal DTR (saturate model)
DTRmax=function(Datain){
       fit2=lm(Y~A1+R+A2+A1*R*A2,data=Datain)        #estimate stage-2 Q-function 
       B=coef(fit2)                                  #theta2_hat
  
       Datain$EQ2max=rep(NA,nrow(Datain))            #create pseudo final outcome for Q1
       for (k in 1:nrow(Datain)){                        
            Datain$EQ2max[k]=Q2max(GA1=Datain$A1[k],GR=Datain$R[k],beta=B)
            }

       fit1=lm(EQ2max~A1,data=Datain)                #estiamte stage-1 Q-function
       G=coef(fit1)                                  #theta1_hat
  
      #estimate optimal DTR
       VQ10=Q1(a1=0,gamma=G)                         #Q1 given A1=0
       VQ11=Q1(a1=1,gamma=G)                         #Q1 given A1=1
       if (VQ10>VQ11) {d1=0} else {d1=1}             #optimal d1
  
       VQ20A=Q2(a1=d1,r=0,a2=0,beta=B)
       VQ21A=Q2(a1=d1,r=0,a2=1,beta=B)
       if (VQ20A>VQ21A) {d20=0} else {d20=1}         #optimal d2 given A1=d1 and R=0
  
       VQ20B=Q2(a1=d1,r=1,a2=0,beta=B)
       VQ21B=Q2(a1=d1,r=1,a2=1,beta=B)
       if (VQ20B>VQ21B) {d21=0} else {d21=1}         #optimal d2 given A1=d1 and R=1
  
      #estimate the value of optimal DTR
       phi0=Q2(a1=d1,r=0,a2=d20,beta=B)              #sequence-specific mean of (d1,R=0,d20)
       phi1=Q2(a1=d1,r=1,a2=d21,beta=B)              #sequence-specific mean of (d1,R=1,d21)
       if (is.nan(phi0)) {phi0=0}
       if (is.nan(phi1)) {phi1=0}
       subD=Datain[which(Datain$A1==d1),]     
       P=table(subD$R)/nrow(subD)                    #P(R=0|d1) and P(R=1|d1)
       value=phi0*P[1]+phi1*P[2]                     #value of optimal DTR E[Y|(d1;d20,d21)]
       return(list(DTR=c(d1,d20,d21),VDTR=as.numeric(value)))
}
DTRmax(Datain=HUB)


#R codes for generating figure 3
#EPI1: stage-1 empirical randomization probabilities (e.g., epi10 epi11) 
EPI1=function(Q10,Q11,sigma1,b){
     Erho10=exp((Q10/sigma1)*log(b))                      #rho_hat
     Erho11=exp((Q11/sigma1)*log(b))
     Epi10=Erho10/(Erho10+Erho11)                         #pi_hat(a1=0|historical data)
     Epi11=Erho11/(Erho10+Erho11)                         #pi_hat(a1=1|historical data)
     return(round(c(Epi10,Epi11),2))
}

#EPI2: stage-2 empirical randomization probabilities (e.g., epi20 epi21|a1,r) 
EPI2=function(Q20,Q21,sigma2,b){
     Erho20=exp((Q20/sigma2)*log(b))                      #rho hat
     Erho21=exp((Q21/sigma2)*log(b))
     Epi20=Erho20/(Erho20+Erho21)                         #pi hat
     Epi21=Erho21/(Erho20+Erho21)
     return(round(c(Epi20,Epi21),2))
}

#API1:Adaptive randomization probabilities: n=n(i)
ARPI1=function(HistPI1=c(0.5,0.5),Q10,Q11,sigma1,b=10,tau=0.75,Nmin,n){
      Erho10=exp((Q10/sigma1)*log(b))                                      #rho_hat
      Erho11=exp((Q11/sigma1)*log(b))
      Epi10=Erho10/(Erho10+Erho11)                                         #pi_hat
      Epi11=Erho11/(Erho10+Erho11)                     
  
      lambda=tau^(1/(b-1))*Nmin/n
      Wrho10=exp(lambda^(b-1)*log(HistPI1[1])+(1-lambda^(b-1))*log(Epi10))    #rho tilde
      Wrho11=exp(lambda^(b-1)*log(HistPI1[2])+(1-lambda^(b-1))*log(Epi11))    
      Wpi10=Wrho10/(Wrho10+Wrho11)                                         #pi tilde
      Wpi11=Wrho11/(Wrho10+Wrho11)              
      return(c(Wpi10,Wpi11))
}

#pih_2: historical pi2
ARPI2=function(HistPI2=c(0.5,0.5),Q20,Q21,sigma2,b=10,tau=0.75,Nmin,n){
      Erho20=exp((Q20/sigma2)*log(b))                                       #rho_hat
      Erho21=exp((Q21/sigma2)*log(b))
      Epi20=Erho20/(Erho20+Erho21)                                          #pi_hat
      Epi21=Erho21/(Erho20+Erho21)
  
      lambda=tau^(1/(b-1))*Nmin/n
      Wrho20=exp(lambda^(b-1)*log(HistPI2[1])+(1-lambda^(b-1))*log(Epi20))    #rho_tilde
      Wrho21=exp(lambda^(b-1)*log(HistPI2[2])+(1-lambda^(b-1))*log(Epi21))
      Wpi20=Wrho20/(Wrho20+Wrho21)                                             #pi_tilde
      Wpi21=Wrho21/(Wrho20+Wrho21)              
      return(round(c(Wpi20,Wpi21),2))
}

#calculate historical randomization probability
histDat=HUB
fit2=lm(Y~A1+R+A2+A1*R*A2,data=histDat)                                 #estimate stage-2 Q-function 
B=coef(fit2)                                                            #theta2_hat
B[which(is.na(B))]=0
SE2=summary(fit2)$sigma                                                 #sigma2_hat

for (k in 1:nrow(histDat)){                                             #create pseudo final outcome for Q1
     histDat$EQ2max[k]=Q2max(GA1=histDat$A1[k],GR=histDat$R[k],beta=B)
}

fit1=lm(EQ2max~A1,data=histDat)                                         #estiamte stage-1 Q-function
G=coef(fit1)                                                            #theta1_hat
SE1=summary(fit1)$sigma                                                 #sigma t=1

VQ10=Q1(a1=0,g=G)
VQ11=Q1(a1=1,g=G)
HPI1=EPI1(Q10=VQ10,Q11=VQ11,sigma1=SE1,b=2)

VQ20A=Q2(a1=0,r=0,a2=0,b=B)                                                #Q2 for A2=0 given a1=0 and r=0
VQ21A=Q2(a1=0,r=0,a2=1,b=B)                                                #Q2 for A2=1 given a1=0 and r=0
HPI2=EPI2(Q20=VQ20A,Q21=VQ21A,sigma2=SE2,b=2)

VQ20B=Q2(a1=0,r=1,a2=0,b=B)                                                #Q2 for A2=0 given a1=0 and r=1
VQ21B=Q2(a1=0,r=1,a2=1,b=B)                                                #Q2 for A2=1 given a1=0 and r=1
HPI2=c(HPI2,EPI2(Q20=VQ20B,Q21=VQ21B,sigma2=SE2,b=2))

VQ20C=Q2(a1=1,r=0,a2=0,b=B)                                                #Q2 for A2=0 given a1=1 and r=0
VQ21C=Q2(a1=1,r=0,a2=1,b=B)                                                #Q2 for A2=1 given a1=1 and r=0
HPI2=c(HPI2,EPI2(Q20=VQ20C,Q21=VQ21C,sigma2=SE2,b=2))

VQ20D=Q2(a1=1,r=1,a2=0,b=B)                                                #Q2 for A2=0 given a1=1 and r=1
VQ21D=Q2(a1=1,r=1,a2=1,b=B)                                                #Q2 for A2=1 given a1=1 and r=1
HPI2=c(HPI2,EPI2(Q20=VQ20D,Q21=VQ21D,sigma2=SE2,b=2))

#simulate sequence-specific outcome of Y
gen.Y=function(Data,mu=c(1,2,3,4,5,6,7,8),sd=3,SEED=10){
      Data$EY[which(Data$A1==0 & Data$R==0 & Data$A2==0)]=mu[1]    
      Data$EY[which(Data$A1==0 & Data$R==0 & Data$A2==1)]=mu[2] 
      Data$EY[which(Data$A1==0 & Data$R==1 & Data$A2==0)]=mu[3]
      Data$EY[which(Data$A1==0 & Data$R==1 & Data$A2==1)]=mu[4]
      Data$EY[which(Data$A1==1 & Data$R==0 & Data$A2==0)]=mu[5]
      Data$EY[which(Data$A1==1 & Data$R==0 & Data$A2==1)]=mu[6]
      Data$EY[which(Data$A1==1 & Data$R==1 & Data$A2==0)]=mu[7]
      Data$EY[which(Data$A1==1 & Data$R==1 & Data$A2==1)]=mu[8]
  
      Data$SDY=sd
      set.seed(SEED)
      Data$Y=rnorm(n=nrow(Data),mean=Data$EY,sd=Data$SDY)
      return(Data)
}

#Simulate a SMART-AR trial (A1,R,A2,Y)
N=400
Nmin=100
Pr=c(0.45,0.5)

seqEY=c(0.86,10.53,11.53,6.04,7.77,4.25,21.60,11.73)

P0=P1=P000=P001=P010=P011=P100=P101=P110=P111=rep(NA,N)
ARPmat1=data.frame(P0,P1)                                                  #output stage-1 AR prob.
ARPmat2=data.frame(P000,P001,P010,P011,P100,P101,P110,P111)                #output stage-2 AR prob.

A1=R=A2=Y=EY=SDY=rep(NA,N)
index=seq(1,N,by=1)
D=data.frame(index,A1,R,A2,EY,SDY,Y)


Pi1=c(0.3,0.7)
Pi2=rep(0.5,8)
for (i in 1:Nmin){
     D$A1[i]=rbinom(n=1, size=1, prob=Pi1)
     if (D$A1[i]==0) {D$R[i]=rbinom(n=1,size=1,prob=Pr[1])} 
               else  {D$R[i]=rbinom(n=1,size=1,prob=Pr[2])} 
     
     if (D$A1[i]==0 & D$R[i]==0) {D$A2[i]=rbinom(n=1,size=1,prob=Pi2[1])} else
     if (D$A1[i]==0 & D$R[i]==1) {D$A2[i]=rbinom(n=1,size=1,prob=Pi2[2])} else
     if (D$A1[i]==1 & D$R[i]==0) {D$A2[i]=rbinom(n=1,size=1,prob=Pi2[3])} else
     if (D$A1[i]==1 & D$R[i]==1) {D$A2[i]=rbinom(n=1,size=1,prob=Pi2[4])}
  
     ARPmat1[i,]=Pi1
     ARPmat2[i,]=Pi2
}
D[1:Nmin,]=gen.Y(Data=D[1:Nmin,],mu=seqEY,sd=5,SEED=12)          #set seed or the first Nmin case
  
#start updating AR probability
j=Nmin+1
while(j<=N){
      print(j)
     #Q-learning to estimate Q1(A1),Q2(A2|A1,R)
      DAT=D[which(D$index<j),]                                   #data for adaptive information
      DAT$EQ2max=NA
      fit2=lm(Y~A1+R+A2+A1*R*A2,data=DAT)                        #estimate stage-2 Q-function 
      B=coef(fit2)                                               #theta2_hat
      SE2=summary(fit2)$sigma                                    #sigma2_hat
  
      for (k in 1:nrow(DAT)){                                    #create pseudo final outcome for Q1
           DAT$EQ2max[k]=Q2max(GA1=DAT$A1[k],GR=DAT$R[k],beta=B)
           }
  
      fit1=lm(EQ2max~A1,data=DAT)                   #estiamte stage-1 Q-function
      G=coef(fit1)                                  #theta1_hat
      SE1=summary(fit1)$sigma                       #sigma t=1
  
     #updated randomization probability for patient j (Pi1,Pi2)
      VQ10=Q1(a1=0,g=G)
      VQ11=Q1(a1=1,g=G)
  
      ni=nrow(DAT)
      newPI1=ARPI1(HistPI1=Pi1,Q10=VQ10,Q11=VQ11,sigma1=SE1,b=2,tau=0.9,Nmin=Nmin,n=ni)
      ARPmat1[j,]=newPI1
  
      VQ20A=Q2(a1=0,r=0,a2=0,b=B)     #Q2 for A2=0 given a1=0 and r=0
      VQ21A=Q2(a1=0,r=0,a2=1,b=B)     #Q2 for A2=1 given a1=0 and r=0
      newPI2=ARPI2(HistPI2=Pi2[1:2],Q20=VQ20A,Q21=VQ21A,sigma2=SE2,b=2,tau=0.9,Nmin=Nmin,n=ni)
  
      VQ20B=Q2(a1=0,r=1,a2=0,b=B)     #Q2 for A2=0 given a1=0 and r=1
      VQ21B=Q2(a1=0,r=1,a2=1,b=B)     #Q2 for A2=1 given a1=0 and r=1
      newPI2=c(newPI2,ARPI2(HistPI2=Pi2[3:4],Q20=VQ20B,Q21=VQ21B,sigma2=SE2,b=2,tau=0.9,Nmin=Nmin,n=ni))
  
      VQ20C=Q2(a1=1,r=0,a2=0,b=B)     #Q2 for A2=0 given a1=1 and r=0
      VQ21C=Q2(a1=1,r=0,a2=1,b=B)     #Q2 for A2=1 given a1=1 and r=0
      newPI2=c(newPI2,ARPI2(HistPI2=Pi2[5:6],Q20=VQ20C,Q21=VQ21C,sigma2=SE2,b=2,tau=0.9,Nmin=Nmin,n=ni))
  
      VQ20D=Q2(a1=1,r=1,a2=0,b=B)     #Q2 for A2=0 given a1=0 and r=0
      VQ21D=Q2(a1=1,r=1,a2=1,b=B)     #Q2 for A2=1 given a1=0 and r=0
      newPI2=c(newPI2,ARPI2(HistPI2=Pi2[7:8],Q20=VQ20D,Q21=VQ21D,sigma2=SE2,b=2,tau=0.9,Nmin=Nmin,n=ni))
      ARPmat2[j,]=newPI2
  
      D$A1[j]=rbinom(n=1,size=1,prob=newPI1[2])
      if (D$A1[j]==0) {D$R[j]=rbinom(n=1,size=1,prob=Pr[1])} 
                else  {D$R[j]=rbinom(n=1,size=1,prob=Pr[2])} 
      if (D$A1[j]==0 & D$R[j]==0) {D$A2[j]=rbinom(n=1,size=1,prob=newPI2[2])} else
      if (D$A1[j]==0 & D$R[j]==1) {D$A2[j]=rbinom(n=1,size=1,prob=newPI2[4])} else
      if (D$A1[j]==1 & D$R[j]==0) {D$A2[j]=rbinom(n=1,size=1,prob=newPI2[6])} else
      if (D$A1[j]==1 & D$R[j]==1) {D$A2[j]=rbinom(n=1,size=1,prob=newPI2[8])} 
  
      D[1:j,]=gen.Y(Data=D[1:j,],mu=seqEY,sd=5,SEED=1)          #set seed or the first Nmin case
      j=j+1
}

par(mfrow=c(2,3))
plot(x=seq(1:nrow(ARPmat1)),y=ARPmat1[,1],type="l",ylim=c(0,1),lwd=2,
     col="blue",xlab="Patients enrolled",ylab="Randomization Probability",
     main="Stage-1 randomization probability")
lines(x=seq(1:nrow(ARPmat1)),y=ARPmat1[,2],lwd=2,col="red")

plot(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,1],type="l",ylim=c(0,1),lwd=2,
     col="blue",xlab="Patients enrolled",ylab="Randomization Probability",
     main="Stage-2 randomization probability")
lines(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,2],lwd=2,col="red")

plot(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,3],type="l",ylim=c(0,1),lwd=2,
     col="blue",xlab="Patients enrolled",ylab="Randomization Probability",
     main="Stage-2 randomization probability")
lines(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,4],lwd=2,col="red")

plot(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,5],type="l",ylim=c(0,1),lwd=2,
     col="blue",xlab="Patients enrolled",ylab="Randomization Probability",
     main="Stage-2 randomization probability")
lines(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,6],lwd=2,col="red")

plot(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,7],type="l",ylim=c(0,1),lwd=2,
     col="blue",xlab="Patients enrolled",ylab="Randomization Probability",
     main="Stage-2 randomization probability")
lines(x=seq(1:nrow(ARPmat2)),y=ARPmat2[,8],lwd=2,col="red")
par(mfrow=c(1,1))


