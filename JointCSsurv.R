list.of.packages <- c("splines2", "numDeriv","statmod")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(splines2)
library(numDeriv)
library(statmod)

Lik.i<-function(t,gam.d,beta.d,alpha.d,c.d,sigma.d,Lij.vec,Rij.vec,DeltaL.vec,DeltaI.vec,Z.vec,Xij.vec,LRij.vec,ibsMat,max.m,loga){
  
  hietai<-prod(sapply(1:length(Lij.vec), 
                      function(o){ifelse(DeltaL.vec[o]==1, 1-exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t)),
                      ifelse(DeltaI.vec[o]==1, (exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))-
                      exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))),
                      exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))))
                      }))*exp(length(Lij.vec)*(sum(Z.vec*alpha.d)+c.d*t))/(1+exp(sum(Z.vec*alpha.d)+c.d*t))^max.m*
    sqrt(2*pi)^(-1)*sigma.d^(-1)*exp(-t^2/2/(sigma.d^2))
  
  if(loga==FALSE){return(hietai)}else if(loga==TRUE){return(log(hietai))}
}

Lik.hat<-function(par,K,P,Q,dat,max.i,LRij.vec,ibsMat,max.m,sigma.vector){
  
  gam.d  <-exp(par[1:K])
  beta.d <-par[(K+1):(K+P)]
  alpha.d<-par[(K+P+1):(K+P+Q)]
  c.d    <-par[K+P+Q+1]
  sigma.d<-exp(par[K+P+Q+2])
  
  lik.store<-0
  for(i in c(1:max.i)){
    dat.ith<-dat[which(dat[,1]==i),]
    if(is.matrix(dat.ith)==FALSE){
      Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
      DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
      Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
    }else{
      Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
      DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
      Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
    }
    
    mode0<-optim(par = 0,fn = Lik.i, method="Brent", lower=-10,upper=10,
                 gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
                 Lij.vec = Lij.vec, Rij.vec = Rij.vec, DeltaL.vec = DeltaL.vec, DeltaI.vec = DeltaI.vec,
                 Z.vec = Z.vec, Xij.vec = Xij.vec, LRij.vec = LRij.vec,
                 ibsMat = ibsMat,max.m = max.m,
                 loga=FALSE,control=list(fnscale=-1,maxit=1000))
    
    sigma0<-1/sqrt(-as.numeric(hessian(func = Lik.i, x=mode0$par,
                                       gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
                                       Lij.vec = Lij.vec, Rij.vec = Rij.vec, DeltaL.vec = DeltaL.vec, DeltaI.vec = DeltaI.vec,
                                       Z.vec = Z.vec, Xij.vec = Xij.vec, LRij.vec = LRij.vec,
                                       ibsMat = ibsMat,max.m = max.m, loga=TRUE)))
    
    out <- gauss.quad(20,"hermite")
    new.node<- mode0$par + sqrt(2)*sigma.vector[i]*out$nodes
    new.wt  <- out$weights*exp(out$nodes^2)
    
    ELik<-sum(sapply(c(1:length(new.node)),function(o){
      sqrt(2)*sigma.vector[i]*Lik.i(new.node[o],gam.d = gam.d, 
                                    beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
                                    Lij.vec = Lij.vec, Rij.vec = Rij.vec, DeltaL.vec = DeltaL.vec, DeltaI.vec = DeltaI.vec,
                                    Z.vec = Z.vec, Xij.vec = Xij.vec, LRij.vec = LRij.vec,
                                    ibsMat = ibsMat,max.m = max.m,loga=FALSE)*new.wt[o]}))
    
    lik.store<-lik.store+log(ELik)
  }
  return(lik.store)
}

Eheta<-function(t,i,P,Q,dat,gam.d,beta.d,alpha.d,c.d,sigma.d,LRij.vec,ibsMat,max.m){
  
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }

  hietai<-prod(sapply(1:length(Lij.vec), 
                      function(o){ifelse(DeltaL.vec[o]==1, 1-exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t)),
                                         ifelse(DeltaI.vec[o]==1, (exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))-
                                                                     exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))),
                                                exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))))
                      }))*exp(length(Lij.vec)*(sum(Z.vec*alpha.d)+c.d*t))/(1+exp(sum(Z.vec*alpha.d)+c.d*t))^max.m*
    sqrt(2*pi)^(-1)*sigma.d^(-1)*exp(-t^2/2/(sigma.d^2))
  return(hietai)
}

logheta<-function(t,i,P,Q,dat,gam.d,beta.d,alpha.d,c.d,sigma.d,LRij.vec,ibsMat,max.m){
  
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  hietai<-prod(sapply(1:length(Lij.vec), 
                      function(o){ifelse(DeltaL.vec[o]==1, 1-exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t)),
                                         ifelse(DeltaI.vec[o]==1, (exp(-sum(gam.d*ibsMat[which(Lij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))-
                                                                     exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))),
                                                exp(-sum(gam.d*ibsMat[which(Rij.vec[o]==LRij.vec),])*exp(sum(Xij.vec[o,]*beta.d)+t))))
                      }))*exp(length(Lij.vec)*(sum(Z.vec*alpha.d)+c.d*t))/(1+exp(sum(Z.vec*alpha.d)+c.d*t))^max.m*
    sqrt(2*pi)^(-1)*sigma.d^(-1)*exp(-t^2/2/(sigma.d^2))
  return(log(hietai))
}

rijk.matrix<-function(i,j,k,K,P,Q,dat,gam.d,beta.d,LRij.vec,ibsMat){
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  lambdaijk<-gam.d[k]*ibsMat[which(Lij.vec[j]==LRij.vec),][k]*exp(sum(Xij.vec[j,]*beta.d))
  Slambdaij <-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j]==LRij.vec),][o]}))*exp(sum(Xij.vec[j,]*beta.d))
  
  return(ifelse(DeltaL.vec[j]==1,lambdaijk/Slambdaij,0))
}

sijk.matrix<-function(i,j,k,K,P,Q,dat,gam.d,beta.d,LRij.vec,ibsMat){
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  epsilonijk<-gam.d[k]*(ibsMat[which(Rij.vec[j]==LRij.vec),][k]-ibsMat[which(Lij.vec[j]==LRij.vec),][k])*exp(sum(Xij.vec[j,]*beta.d))
  Sepsilonij<-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j]==LRij.vec),][o]-ibsMat[which(Lij.vec[j]==LRij.vec),][o])}))*exp(sum(Xij.vec[j,]*beta.d))
  
  return(ifelse(DeltaI.vec[j]==1,epsilonijk/Sepsilonij,0))
}

E.Y.matrix<-function(i,j,k,K,P,Q,dat,gam.d,beta.d,alpha.d,c.d,LRij.vec,ibsMat,sigma.vector,in.pow,out.pow,a,b,c,ij,new.node,new.wt,HIETAI){
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  lambdaijk<-gam.d[k]*ibsMat[which(Lij.vec[j]==LRij.vec),][k]*exp(sum(Xij.vec[j,]*beta.d))
  Slambdaij<-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j]==LRij.vec),][o]}))*exp(sum(Xij.vec[j,]*beta.d))
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  
  if(in.pow==1){
    if(ij==TRUE){Y<-pmax(0,ifelse(rep(DeltaL.vec[j]==1,length(new.node)),(exp(new.node)*Slambdaij/(1-exp(-exp(new.node)*Slambdaij)))^out.pow,0))
    }else if(ij==FALSE){Y<-pmax(0,ifelse(rep(DeltaL.vec[j]==1,length(new.node)),(exp(new.node)*lambdaijk/(1-exp(-exp(new.node)*Slambdaij)))^out.pow,0))}
    EY.d<-sqrt(2)*sigma.vector[i]*(Y*(new.node^a)*exp(b*new.node)*(1+exp(sum(-Z.vec*alpha.d)-c.d*new.node))^(-c)*HIETAI[i,])%*%new.wt
  }else if(in.pow==2){
    Y<-pmax(0,ifelse(rep(DeltaL.vec[j]==1,length(new.node)),((exp(new.node)*Slambdaij)^2+exp(new.node)*Slambdaij)/(1-exp(-exp(new.node)*Slambdaij)),0))
    EY.d<-sqrt(2)*sigma.vector[i]*(Y*HIETAI[i,])%*%new.wt
  }
  return(EY.d/Eh.d)
}

E.W.matrix<-function(i,j,k,K,P,Q,dat,gam.d,beta.d,alpha.d,c.d,LRij.vec,ibsMat,sigma.vector,in.pow,out.pow,a,b,c,ij,new.node,new.wt,HIETAI){
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  epsilonijk<-gam.d[k]*(ibsMat[which(Rij.vec[j]==LRij.vec),][k]-ibsMat[which(Lij.vec[j]==LRij.vec),][k])*exp(sum(Xij.vec[j,]*beta.d))
  Sepsilonij <-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j]==LRij.vec),][o]-ibsMat[which(Lij.vec[j]==LRij.vec),][o])}))*exp(sum(Xij.vec[j,]*beta.d))
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  
  if(in.pow==1){
    if(ij==TRUE){W<-pmax(0,ifelse(rep(DeltaI.vec[j]==1,length(new.node)),(exp(new.node)*Sepsilonij/(1-exp(-exp(new.node)*Sepsilonij)))^out.pow,0))
    }else if(ij==FALSE){W<-pmax(0,ifelse(rep(DeltaI.vec[j]==1,length(new.node)),(exp(new.node)*epsilonijk/(1-exp(-exp(new.node)*Sepsilonij)))^out.pow,0))}
    EW.d<-sqrt(2)*sigma.vector[i]*(W*(new.node^a)*exp(b*new.node)*(1+exp(sum(-Z.vec*alpha.d)-c.d*new.node))^(-c)*HIETAI[i,])%*%new.wt
  }else if(in.pow==2){
    W<-pmax(0,ifelse(rep(DeltaI.vec[j]==1,length(new.node)),((exp(new.node)*Sepsilonij)^2+exp(new.node)*Sepsilonij)/(1-exp(-exp(new.node)*Sepsilonij)),0))
    EW.d<-sqrt(2)*sigma.vector[i]*(W*HIETAI[i,])%*%new.wt
  }
  return(EW.d/Eh.d)
}

E.YW.matrix<-function(i,j1,j2,k1,k2,K,P,Q,dat,gam.d,beta.d,LRij.vec,ibsMat,sigma.vector,a,b,ij,new.node,new.wt,HIETAI){
  
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  lambdaijk1<-gam.d[k1]*ibsMat[which(Lij.vec[j1]==LRij.vec),][k1]*exp(sum(Xij.vec[j1,]*beta.d))
  lambdaijk2<-gam.d[k2]*ibsMat[which(Lij.vec[j2]==LRij.vec),][k2]*exp(sum(Xij.vec[j2,]*beta.d))
  epsilonijk1<-gam.d[k1]*(ibsMat[which(Rij.vec[j1]==LRij.vec),][k1]-ibsMat[which(Lij.vec[j1]==LRij.vec),][k1])*exp(sum(Xij.vec[j1,]*beta.d))
  epsilonijk2<-gam.d[k2]*(ibsMat[which(Rij.vec[j2]==LRij.vec),][k2]-ibsMat[which(Lij.vec[j2]==LRij.vec),][k2])*exp(sum(Xij.vec[j2,]*beta.d))
  Slambdaij1 <-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j1]==LRij.vec),][o]}))*exp(sum(Xij.vec[j1,]*beta.d))
  Slambdaij2 <-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j2]==LRij.vec),][o]}))*exp(sum(Xij.vec[j2,]*beta.d))
  Sepsilonij1 <-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j1]==LRij.vec),][o]-ibsMat[which(Lij.vec[j1]==LRij.vec),][o])}))*exp(sum(Xij.vec[j1,]*beta.d))
  Sepsilonij2 <-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j2]==LRij.vec),][o]-ibsMat[which(Lij.vec[j2]==LRij.vec),][o])}))*exp(sum(Xij.vec[j2,]*beta.d))
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  if(ij==TRUE){
    Y1<-pmax(0,ifelse(rep(DeltaL.vec[j1]==1,length(new.node)),exp(new.node)*Slambdaij1/(1-exp(-exp(new.node)*Slambdaij1)),0))
    W1<-pmax(0,ifelse(rep(DeltaI.vec[j1]==1,length(new.node)),exp(new.node)*Sepsilonij1/(1-exp(-exp(new.node)*Sepsilonij1)),0))
    Y2<-pmax(0,ifelse(rep(DeltaL.vec[j2]==1,length(new.node)),exp(new.node)*Slambdaij2/(1-exp(-exp(new.node)*Slambdaij2)),0))
    W2<-pmax(0,ifelse(rep(DeltaI.vec[j2]==1,length(new.node)),exp(new.node)*Sepsilonij2/(1-exp(-exp(new.node)*Sepsilonij2)),0))
  }else if(ij==FALSE){
    Y1<-pmax(0,ifelse(rep(DeltaL.vec[j1]==1,length(new.node)),exp(new.node)*lambdaijk1/(1-exp(-exp(new.node)*Slambdaij1)),0))
    W1<-pmax(0,ifelse(rep(DeltaI.vec[j1]==1,length(new.node)),exp(new.node)*epsilonijk1/(1-exp(-exp(new.node)*Sepsilonij1)),0))
    Y2<-pmax(0,ifelse(rep(DeltaL.vec[j2]==1,length(new.node)),exp(new.node)*lambdaijk2/(1-exp(-exp(new.node)*Slambdaij2)),0))
    W2<-pmax(0,ifelse(rep(DeltaI.vec[j2]==1,length(new.node)),exp(new.node)*epsilonijk2/(1-exp(-exp(new.node)*Sepsilonij2)),0))}
  
  if(a==2){EY.d<-sqrt(2)*sigma.vector[i]*(Y1*Y2*HIETAI[i,])%*%new.wt
  }else if(b==2){EY.d<-sqrt(2)*sigma.vector[i]*(W1*W2*HIETAI[i,])%*%new.wt
  }else{EY.d<-sqrt(2)*sigma.vector[i]*(Y1*W2*HIETAI[i,])%*%new.wt}
  
  return(EY.d/Eh.d)
}

E.YW.sp.matrix<-function(i,j1,j2,k1,k2,K,P,Q,dat,gam.d,beta.d,LRij.vec,ibsMat,sigma.vector,a,b,YW,new.node,new.wt,HIETAI){
  
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  Slambdaij1 <-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j1]==LRij.vec),][o]}))*exp(sum(Xij.vec[j1,]*beta.d))
  Sepsilonij1 <-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j1]==LRij.vec),][o]-ibsMat[which(Lij.vec[j1]==LRij.vec),][o])}))*exp(sum(Xij.vec[j1,]*beta.d))
  lambdaijk2 <-gam.d[k2]*ibsMat[which(Lij.vec[j2]==LRij.vec),][k2]*exp(sum(Xij.vec[j2,]*beta.d))
  epsilonijk2<-gam.d[k2]*(ibsMat[which(Rij.vec[j2]==LRij.vec),][k2]-ibsMat[which(Lij.vec[j2]==LRij.vec),][k2])*exp(sum(Xij.vec[j2,]*beta.d))
  Slambdaij2 <-sum(sapply(1:K, function(o){gam.d[o]*ibsMat[which(Lij.vec[j2]==LRij.vec),][o]}))*exp(sum(Xij.vec[j2,]*beta.d))
  Sepsilonij2 <-sum(sapply(1:K, function(o){gam.d[o]*(ibsMat[which(Rij.vec[j2]==LRij.vec),][o]-ibsMat[which(Lij.vec[j2]==LRij.vec),][o])}))*exp(sum(Xij.vec[j2,]*beta.d))
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  
  Y1<-pmax(0,ifelse(rep(DeltaL.vec[j1]==1,length(new.node)),exp(new.node)*Slambdaij1/(1-exp(-exp(new.node)*Slambdaij1)),0))
  W1<-pmax(0,ifelse(rep(DeltaI.vec[j1]==1,length(new.node)),exp(new.node)*Sepsilonij1/(1-exp(-exp(new.node)*Sepsilonij1)),0))
  Y2<-pmax(0,ifelse(rep(DeltaL.vec[j2]==1,length(new.node)),exp(new.node)*lambdaijk2/(1-exp(-exp(new.node)*Slambdaij2)),0))
  W2<-pmax(0,ifelse(rep(DeltaI.vec[j2]==1,length(new.node)),exp(new.node)*epsilonijk2/(1-exp(-exp(new.node)*Sepsilonij2)),0))
  
  if(a==2){EY.d<-sqrt(2)*sigma.vector[i]*(Y1*Y2*HIETAI[i,])%*%new.wt
  }else if(b==2){EY.d<-sqrt(2)*sigma.vector[i]*(W1*W2*HIETAI[i,])%*%new.wt
  }else if((a==1)&(b==1)&(YW==TRUE)){EY.d<-sqrt(2)*sigma.vector[i]*(Y1*W2*HIETAI[i,])%*%new.wt
  }else if((a==1)&(b==1)&(YW==FALSE)){EY.d<-sqrt(2)*sigma.vector[i]*(W1*Y2*HIETAI[i,])%*%new.wt}
  
  return(EY.d/Eh.d)
}

etai.abcd.vector<-function(i,P,Q,dat,alpha.d,c.d,sigma.vector,a,b,c,d,new.node,new.wt,HIETAI){
  
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  Eetai.i<-sqrt(2)*sigma.vector[i]*((new.node^a)*exp(b*new.node)*(1+exp(sum(-Z.vec*alpha.d)-c.d*new.node))^(-c)*
                                      exp(sum(-Z.vec*alpha.d)-c.d*new.node)^(d)*HIETAI[i,])%*%new.wt
  
  return(Eetai.i/Eh.d)
}

mode.value.i<-function(i,P,Q,dat,gam.d,beta.d,alpha.d,c.d,sigma.d,LRij.vec,ibsMat,max.m){
  
  mode0<-optim(par = 0,fn = Eheta, method="Brent", lower=-10,upper=10, i=i, P=P, Q=Q,
               dat = dat,gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
               ibsMat = ibsMat, max.m = max.m, LRij.vec=LRij.vec,
               control=list(fnscale=-1,maxit=1000))
  
  return(mode0$par)
}

sigma.value.i<-function(i,P,Q,dat,gam.d,beta.d,alpha.d,c.d,sigma.d,LRij.vec,ibsMat,max.m){
  
  mode0<-optim(par = 0,fn = Eheta, method="Brent", lower=-10,upper=10, i=i, P=P, Q=Q,
               dat=dat,gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
               ibsMat = ibsMat, max.m = max.m, LRij.vec=LRij.vec,
               control=list(fnscale=-1,maxit=1000))
  
  sigma0<-1/sqrt(-as.numeric(hessian(func = logheta, x=mode0$par, i=i, P=P, Q=Q,
                                     dat=dat,gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
                                     ibsMat = ibsMat, max.m = max.m, LRij.vec=LRij.vec)))
  
  return(sigma0)
}

gam.k<-function(beta,k,P,Q,dat,T010,Yijk000,Wijk000,LRij.vec,max.i,ibsMat){
  store0<-0
  store1<-0
  for(i in c(1:max.i)){
    dat.ith<-dat[which(dat[,1]==i),]
    if(is.matrix(dat.ith)==FALSE){
      Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
      DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
      Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
    }else{
      Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
      DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
      Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
    }
    
    for(q in c(1:length(Lij.vec))){
      store0<-store0+(Yijk000[q,k,i]+Wijk000[q,k,i])
      store1<-store1+exp(sum(beta*Xij.vec[q,]))*T010[i]*ifelse(DeltaL.vec[q]==1,ibsMat[which(Lij.vec[q]==LRij.vec),][k], ibsMat[which(Rij.vec[q]==LRij.vec),][k])
    }}
  return(store0/store1)
}

Elog.ac<-function(i,P,Q,dat,T100,max.m,alpha,c,sigma.vector,new.node,new.wt,HIETAI){
  dat.ith<-dat[which(dat[,1]==i),]
  if(is.matrix(dat.ith)==FALSE){
    Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
    DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
    Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
  }else{
    Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
    DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
    if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
    Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
  }
  
  Eh.d<-sqrt(2)*sigma.vector[i]*(HIETAI[i,]%*%new.wt)
  Elog.ac.d<-sqrt(2)*sigma.vector[i]*(log(1+exp(sum(Z.vec*alpha)+c*new.node))*HIETAI[i,])%*%new.wt
  
  return(length(Lij.vec)*(sum(Z.vec*alpha)+c*T100[i])-max.m*Elog.ac.d/Eh.d)
}

M2<-function(par,P,Q,dat,T100,max.i,max.m,sigma.vector,node.mat,new.wt,HIETAI){
  sum(sapply(c(1:max.i),function(o){Elog.ac(i = o, dat=dat, P=P, Q=Q,max.m = max.m,alpha = par[1:(1+Q-1)],c = par[(1+Q)],
                                            new.node = node.mat[o,],T100=T100,
                                            sigma.vector=sigma.vector, new.wt = new.wt,HIETAI = HIETAI)}))
}

JointCSsurvSIM<-function(seed = NA, n, m, beta, alpha, kappa, sigma){
  if(is.numeric(seed)){set.seed(seed)}
  #cluster-specific variable
  Z<-eta<-cs<-NULL
  repeat{
    A   <-rnorm(n = 1, mean = 0, sd = 1)      #Z_i1
    B   <-rnorm(n = 1, mean = 0, sd = sigma)  #eta_i
    C   <-rbinom(n = 1, size = m, prob = exp(alpha[1]+A*alpha[2]+kappa*B)/(1+exp(alpha[1]+A*alpha[2]+kappa*B))) #N_i
    if((C>0)&(C<=m)){Z<-c(Z,A); eta<-c(eta,B); cs<-c(cs,C)}
    if(length(Z)==n){break}
  }
  
  #interval-censored failure time
  store<-NULL
  for(i in c(1:n)){
    X        <-rnorm(n = cs[i], mean = 0, sd = 1)        
    Sij      <-runif(n = cs[i],min = 0,max = 1)
    regressor<-exp(X*beta+eta[i])
    tij      <-((-log(Sij)/regressor)/0.25)^(1/2)
    ct       <- 4
    
    checkup.times<-lapply(1:cs[i],function(w){x<-cumsum(runif(51,min = 0.1,max = 1))
    x<-c(0,x[which(x<ct+2)])
    return(x)})
    
    which.interval<-sapply(1:cs[i],function(x){findInterval(tij[x],checkup.times[[x]])})
    
    Lij1<-ifelse(tij>=ct, NA, sapply(1:cs[i],function(y){checkup.times[[y]][which.interval[y]]}))
    Rij<-ifelse(tij>=ct, ct, pmin(ct,sapply(1:cs[i],function(z){checkup.times[[z]][which.interval[z]+1]})))
    
    Lij<-ifelse((!is.na(Lij1))&(Lij1==0), Rij, Lij1)
    Rij<-ifelse((!is.na(Lij))&(!is.na(Rij))&(Lij==Rij),NA,Rij)
    
    DeltaL<-ifelse(is.na(Rij),1,0)
    DeltaI<-ifelse(!is.na(Lij)&!is.na(Rij),1,0)
    
    #Data structure by column: identifier, cluster size, left and right endpoints of an interval, 
    #                          censoring indicators for left- and interval-censored, covariates
    #Remarks: We coded (Lij,NA) and (NA,Rij) for left- and right-censored observations, respectively
    store<-rbind(store,cbind(i,cs[i],Lij,Rij,DeltaL,DeltaI,X,Z[i]))
  }
  
  store<-data.frame(store)
  names(store)<-c("id","cs","Lij","Rij","DL","DI","X","Z")
  
  return(store)
}


#gam_0=NA; beta_0=NA; alpha_0=NA; kappa_0=NA; sigma_0=NA; TRACE=TRUE

JointCSsurvEST<-function(data, K=7, P, Q, deg=3, max.m, M=20, tolerance=10^{-3}, 
                         gam_0=NA, beta_0=NA, alpha_0=NA, kappa_0=NA, sigma_0=NA, TRACE=FALSE){

  part   <-K-deg+1
  LRij  <- unique(c(data$Lij,data$Rij))
  LRij  <- c(0,sort(unique(LRij[(LRij>=0)&(LRij<=max(LRij,na.rm=TRUE))])))
  if(K<deg){print("ERROR: K cannot be less than the degree of polynomial splines (deg).");break}
  if(K==deg){knots<-NULL}else{knots <- quantile(LRij,probs = seq(1/part,(part-1)/part,1/part),na.rm = T)}  #adaptive
  ibsMat<- iSpline(LRij, knots = knots, degree = deg, intercept = FALSE)
  
  max.i<-length(unique(data$id))
  
  if(is.numeric(gam_0)==FALSE){gam.d <- rep(2,K)}else
    if(length(gam_0)!=K){print("ERROR: The length of gam_0 does not equal to K.");break}else{gam.d <- gam_0}
  if(is.numeric(beta_0)==FALSE){beta.d <- rep(0,P)}else
    if(length(beta_0)!=P){print("ERROR: The length of beta_0 does not equal to P.");break}else{beta.d <- beta_0}
  if(is.numeric(alpha_0)==FALSE){alpha.d <- rep(0,Q)}else 
    if(length(alpha_0)!=Q){print("ERROR: The length of alpha_0 does not equal to Q.");break}else{alpha.d <- alpha_0}
  if(is.numeric(kappa_0)==FALSE){c.d <- 0}else 
    if(length(kappa_0)!=1){print("ERROR: kappa_0 is not a constant.");break}else{c.d <- kappa_0}
  if(is.numeric(sigma_0)==FALSE){sigma.d <- 2}else
    if((length(sigma_0)!=1)|(sigma_0<=0)){print("ERROR: sigma_0 is not a positive constant.");break}else{sigma.d <- sigma_0}

  I<-0
  tick<-0
  repeat{
    I<-I+1
    if(tick==0){theta0<-c(gam.d,beta.d,alpha.d,c.d,sigma.d)}
    if(tick==2){
      r<-(theta1-theta0)
      v<-(theta2-theta1)-r
      a<- -sqrt(sum(r^2))/sqrt(sum(v^2))
      theta.dash<-theta0-2*a*r+a^2*v
      
      if((prod(theta.dash[c(1:K,K+P+Q+2)]>0))&(theta.dash[K+P+Q+2]>0.1)){
        gam.d  <-theta.dash[1:K]            ; beta.d<-theta.dash[(K+1):(K+P)];
        alpha.d<-theta.dash[(K+P+1):(K+P+Q)]; c.d   <-theta.dash[(K+P+Q+1)];
        sigma.d<-theta.dash[(K+P+Q+2)]
      }}
    
    mode.vector<-sapply(c(1:max.i),function(o){mode.value.i(i=o, P=P, Q=Q, dat=as.matrix(data), LRij.vec=LRij, ibsMat = ibsMat, max.m=max.m,
                                                            gam.d=gam.d, beta.d=beta.d, alpha.d=alpha.d, c.d=c.d, sigma.d=sigma.d)})
    
    sigma.vector<-sapply(c(1:max.i),function(o){sigma.value.i(o, P=P, Q=Q, dat=as.matrix(data), LRij.vec=LRij, ibsMat = ibsMat, max.m=max.m,
                                                              gam.d=gam.d, beta.d=beta.d, alpha.d=alpha.d, c.d=c.d, sigma.d=sigma.d)})
    
    beta.dummy<-beta.d
    
    out     <- gauss.quad(M,"hermite")
    new.wt  <- out$weights*exp(out$nodes^2)
    node.mat<- t(sapply(c(1:max.i),function(i){mode.vector[i] + sqrt(2)*sigma.vector[i]*out$nodes}))
    
    HIETAI<-matrix(0, nrow=max.i, ncol=M)
    for(i in c(1:max.i)){
      HIETAI[i,]<-sapply(c(1:M),function(o){Eheta(t = node.mat[i,o],i = i,P = P,Q = Q,dat = as.matrix(data),
                                                  gam.d = gam.d, beta.d = beta.d, alpha.d = alpha.d, c.d = c.d, sigma.d = sigma.d,
                                                  LRij.vec = LRij,ibsMat = ibsMat,max.m = max.m)})
    }
    
    T010<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=1,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    T100<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    T200<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    
    Yij000<-Wij000<-Lambda.Lij<-Lambda.Rij<-array(0,c(max.i,max.m))
    for(i in c(1:max.i)){
      dat.ith<-as.matrix(data)[which(data[,1]==i),]
      if(is.matrix(dat.ith)==FALSE){Lij.vec<-dat.ith[3];Rij.vec<-dat.ith[4]
      }else{Lij.vec<-dat.ith[,3];Rij.vec<-dat.ith[,4]}
      new.node<- mode.vector[i] + sqrt(2)*sigma.vector[i]*out$nodes
      for(j in c(1:length(which(data$id==i)))){
        Yij000[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                                dat=as.matrix(data), gam.d=gam.d,beta.d=beta.d,alpha.d=alpha.d,c.d=c.d,
                                LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wij000[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                                dat=as.matrix(data), gam.d=gam.d,beta.d=beta.d,alpha.d=alpha.d,c.d=c.d,
                                LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Lambda.Lij[i,j]<-sum(gam.d*ibsMat[which(Lij.vec[j]==LRij),])
        Lambda.Rij[i,j]<-sum(gam.d*ibsMat[which(Rij.vec[j]==LRij),])
      }}
    
    score.beta<-sapply(1:P, function(p){
      sum(sapply(c(1:max.i),function(i){
        dat.ith<-as.matrix(data)[which(data[,1]==i),]
        if(is.matrix(dat.ith)==FALSE){
          Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
          DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
          if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
          Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
        }else{
          Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
          DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
          if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
          Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
        }
        store00<-0
        for(j in c(1:length(DeltaL.vec))){
          store00<-store00+Xij.vec[j,p]*(Yij000[i,j]+Wij000[i,j])-
            Xij.vec[j,p]*exp(sum(Xij.vec[j,]*beta.d))*T010[i]*ifelse(DeltaL.vec[j]==1,Lambda.Lij[i,j],Lambda.Rij[i,j])} 
        return(store00)
      })) })
    
    hess.beta<-matrix(0,nrow=P,ncol=P)
    for(p1 in c(1:P)){
      for(p2 in c(1:P)){
        hess.beta[p1 ,p2]<- sum(sapply(c(1:max.i),function(i){
          dat.ith<-as.matrix(data)[which(data[,1]==i),]
          if(is.matrix(dat.ith)==FALSE){
            Lij.vec<-dat.ith[3];   Rij.vec<-dat.ith[4];
            DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
            if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
            Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
          }else{
            Lij.vec<-dat.ith[,3];   Rij.vec<-dat.ith[,4]; 
            DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
            if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
            Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
          }
          store00<-0
          for(j in c(1:length(DeltaL.vec))){
            store00<-store00-Xij.vec[j,p1]*Xij.vec[j,p2]*exp(sum(Xij.vec[j,]*beta.d))*T010[i]*ifelse(DeltaL.vec[j]==1,Lambda.Lij[i,j],Lambda.Rij[i,j])} 
          return(store00)
        }))
      }}
    
    beta.d<-beta.d-as.numeric(solve(hess.beta)%*%score.beta)
    d1    <-beta.d-beta.dummy
    
    gam.dummy<-gam.d
    T010<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=1,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    Yijk000<-Wijk000<-array(0,c(max.m,K,max.i)) #row, #col, #matrix
    for(i in c(1:max.i)){
      new.node<-node.mat[i,]
      for(j in c(1:length(which(data$id==i)))){
        for(k in c(1:K)){
          Yijk000[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=FALSE,
                                     dat=as.matrix(data), gam.d=gam.d,beta.d=beta.d,alpha.d=alpha.d,c.d=c.d,
                                     LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
          Wijk000[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=FALSE,
                                     dat=as.matrix(data), gam.d=gam.d,beta.d=beta.d,alpha.d=alpha.d,c.d=c.d,
                                     LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        }}}
    for(g2 in c(1:K)){
      gam.d[g2]<-gam.k(beta=beta.d,k=g2,P=P,Q=Q,dat=as.matrix(data),T010=T010,Yijk000=Yijk000,Wijk000=Wijk000,LRij.vec=LRij,max.i=max.i,ibsMat=ibsMat)}
    d0<-gam.d-gam.dummy
    
    #Update alpha, c
    T100<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    alpha.c.new<-optim(par = c(alpha.d,c.d),fn = M2,gr = NULL,method = "BFGS",
                       P=P, Q=Q,dat=as.matrix(data),T100=T100, max.i = max.i, max.m = max.m,
                       sigma.vector=sigma.vector,node.mat=node.mat,new.wt=new.wt,HIETAI=HIETAI,
                       control=list(fnscale=-1,maxit=1000,reltol=10^-8))
    
    d2<-alpha.c.new$par[1:Q]-alpha.d
    d3<-alpha.c.new$par[(Q+1)]-c.d
    alpha.d<-alpha.c.new$par[1:Q]
    c.d    <-alpha.c.new$par[(Q+1)]
    
    T200<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                         alpha.d=alpha.d,c.d = c.d,sigma.vector=sigma.vector,HIETAI=HIETAI)})
    #Update sigma
    d4<- sqrt((max.i)^(-1)*sum(T200)) - sigma.d
    sigma.d<-sqrt((max.i)^(-1)*sum(T200))
    
    if(tick==0){theta1<-c(gam.d,beta.d,alpha.d,c.d,sigma.d)}
    if(tick==1){theta2<-c(gam.d,beta.d,alpha.d,c.d,sigma.d)}
    tick<-tick+1
    if(tick==3){tick=0}
    
    logLike<-Lik.hat(par = c(log(gam.d),beta.d,alpha.d,c.d,log(sigma.d)),sigma.vector=sigma.vector,
                     K=K,P=P,Q=Q,dat=as.matrix(data),max.i=max.i,LRij.vec=LRij,ibsMat=ibsMat,max.m=max.m)
    
    if(TRACE==TRUE){print(round(c(I,alpha.d,c.d,beta.d,sigma.d,max(abs(c(d0,d1,d2,d3,d4))),logLike),4))}
    
    if(max(abs(c(d0,d1,d2,d3,d4)))<tolerance){
      para     <-c(log(gam.d),beta.d,alpha.d,c.d,log(sigma.d))
      gam.hat  <-gam.d
      beta.hat <-beta.d
      alpha.hat<-alpha.d
      c.hat    <-c.d
      sigma.hat<-sigma.d
      dist     <-max(abs(c(d0,d1,d2,d3,d4)))
      break
    }
  }
  
  #Interval Estimation
  HIETAI<-matrix(0, nrow=max.i, ncol=M)
  for(i in c(1:max.i)){
    HIETAI[i,]<-sapply(c(1:M),function(o){Eheta(t = node.mat[i,o],i = i,P = P,Q = Q,dat = as.matrix(data),
                                                gam.d = gam.hat, beta.d = beta.hat, alpha.d = alpha.hat, c.d = c.hat, sigma.d = sigma.hat,
                                                LRij.vec = LRij,ibsMat = ibsMat,max.m = max.m)})
  }
  mode.vector<-sapply(c(1:max.i),function(o){mode.value.i(i=o, P=P, Q=Q, dat=as.matrix(data), LRij.vec=LRij, ibsMat = ibsMat, max.m=max.m,
                                                          gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,sigma.d=sigma.hat)})
  sigma.vector<-sapply(c(1:max.i),function(o){sigma.value.i(o, P=P, Q=Q, dat=as.matrix(data), LRij.vec=LRij, ibsMat = ibsMat, max.m=max.m,
                                                            gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,sigma.d=sigma.hat)})
  
  Yij000<-Wij000<-Yij001<-Wij001<-Yij010<-Wij010<-Yij2<-Wij2<-Lambda.Lij<-Lambda.Rij<-EYij2<-EWij2<-array(0,c(max.i,max.m)) #EYij, EWij's dimension (#cluster*max(j))
  Yij100<-Wij100<-Yij101<-Wij101<-Yij200<-Wij200<-array(0,c(max.i,max.m))
  for(i in c(1:max.i)){
    dat.ith<-as.matrix(data)[which(data[,1]==i),]
    if(is.matrix(dat.ith)==FALSE){Lij.vec<-dat.ith[3];Rij.vec<-dat.ith[4]
    }else{Lij.vec<-dat.ith[,3];Rij.vec<-dat.ith[,4]}
    new.node<-node.mat[i,]
    for(j in c(1:length(which(data$id==i)))){
      Yij000[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij000[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij001[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=1,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij001[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=1,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij010[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=1,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij010[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=1,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij100[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij100[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij101[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=1,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij101[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=1,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij200[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=2,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij200[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=2,b=0,c=0,ij=TRUE,
                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                              LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      EYij2[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=2,a=0,b=0,c=0,ij=TRUE,
                             dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                             LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      EWij2[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=1,out.pow=2,a=0,b=0,c=0,ij=TRUE,
                             dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                             LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Yij2[i,j]<-E.Y.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=2,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                            dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                            LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Wij2[i,j]<-E.W.matrix(i=i,j=j,k=1,K=K,P=P,Q=Q,in.pow=2,out.pow=1,a=0,b=0,c=0,ij=TRUE,
                            dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                            LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      Lambda.Lij[i,j]<-sum(gam.hat*ibsMat[which(Lij.vec[j]==LRij),])
      Lambda.Rij[i,j]<-sum(gam.hat*ibsMat[which(Rij.vec[j]==LRij),])
    }}
  
  
  Yijk000<-Wijk000<-Yijk010<-Wijk010<-rijk<-sijk<-Yijk2<-Wijk2<-EYijk2<-EWijk2<-array(0,c(max.m,K,max.i)) #EYijk dimension (#cluster*max(j)*K)
  Yijk001<-Wijk001<-Yijk100<-Wijk100<-Yijk101<-Wijk101<-Yijk200<-Wijk200<-array(0,c(max.m,K,max.i))
  for(i in c(1:max.i)){
    new.node<-node.mat[i,]
    for(j in c(1:length(which(data$id==i)))){
      for(k in c(1:K)){
        Yijk000[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk000[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Yijk001[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=1,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk001[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=0,c=1,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Yijk010[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=1,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk010[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=0,b=1,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Yijk100[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk100[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Yijk101[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=1,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk101[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=1,b=0,c=1,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Yijk200[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=2,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        Wijk200[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=1,a=2,b=0,c=0,ij=FALSE,
                                   dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                   LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        EYijk2[j,k,i]<-E.Y.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=2,a=0,b=0,c=0,ij=FALSE,
                                  dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                  LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        EWijk2[j,k,i]<-E.W.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,in.pow=1,out.pow=2,a=0,b=0,c=0,ij=FALSE,
                                  dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,alpha.d=alpha.hat,c.d=c.hat,
                                  LRij.vec = LRij,ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        rijk[j,k,i]<-rijk.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,ibsMat = ibsMat)
        sijk[j,k,i]<-sijk.matrix(i=i,j=j,k=k,K=K,P=P,Q=Q,dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,ibsMat = ibsMat)
        Yijk2[j,k,i]<-rijk[j,k,i]*(1-rijk[j,k,i])*Yij000[i,j]+rijk[j,k,i]^2*(Yij2[i,j]-EYij2[i,j])+EYijk2[j,k,i]
        Wijk2[j,k,i]<-sijk[j,k,i]*(1-sijk[j,k,i])*Wij000[i,j]+sijk[j,k,i]^2*(Wij2[i,j]-EWij2[i,j])+EWijk2[j,k,i]
      }}}
  
  YikYil<-WikWil<-YikWil<-array(0,c(max.m,max.m,max.i))
  for(i in c(1:max.i)){
    new.node<-node.mat[i,]
    for(j1 in c(1:length(which(data$id==i)))){
      for(j2 in c(1:length(which(data$id==i)))){
        YikYil[j1,j2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=1,k2=1,K=K,P=P,Q=Q,a=2,b=0,ij=TRUE,
                                     dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                     ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        WikWil[j1,j2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=1,k2=1,K=K,P=P,Q=Q,a=0,b=2,ij=TRUE,
                                     dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                     ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        YikWil[j1,j2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=1,k2=1,K=K,P=P,Q=Q,a=1,b=1,ij=TRUE,
                                     dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                     ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
      }}}
  
  YikYijl<-WikWijl<-YikWijl<-WikYijl<-array(0,c(max.m,max.m,K,max.i))
  for(i in c(1:max.i)){
    new.node<-node.mat[i,]
    for(j1 in c(1:length(which(data$id==i)))){
      for(j2 in c(1:length(which(data$id==i)))){
        for(k2 in c(1:K)){
          YikYijl[j1,j2,k2,i]<-E.YW.sp.matrix(i=i,j1=j1,j2=j2,k1=1,k2=k2,K=K,P=P,Q=Q,a=2,b=0,YW=FALSE,
                                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                              ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
          WikWijl[j1,j2,k2,i]<-E.YW.sp.matrix(i=i,j1=j1,j2=j2,k1=1,k2=k2,K=K,P=P,Q=Q,a=0,b=2,YW=FALSE,
                                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                              ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
          YikWijl[j1,j2,k2,i]<-E.YW.sp.matrix(i=i,j1=j1,j2=j2,k1=1,k2=k2,K=K,P=P,Q=Q,a=1,b=1,YW=TRUE,
                                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                              ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
          WikYijl[j1,j2,k2,i]<-E.YW.sp.matrix(i=i,j1=j1,j2=j2,k1=1,k2=k2,K=K,P=P,Q=Q,a=1,b=1,YW=FALSE,
                                              dat=as.matrix(data), gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                              ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
        }}}}
  
  
  YijkYijl<-WijkWijl<-YijkWijl<-array(0,c(max.m,max.m,K,K,max.i))
  for(i in c(1:max.i)){
    new.node<-node.mat[i,]
    for(j1 in c(1:length(which(data$id==i)))){
      for(j2 in c(1:length(which(data$id==i)))){
        for(k1 in c(1:K)){
          for(k2 in c(1:K)){
            YijkYijl[j1,j2,k1,k2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=k1,k2=k2,K=K,P=P,Q=Q,a=2,b=0,ij=FALSE,
                                                 dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                                 ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
            WijkWijl[j1,j2,k1,k2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=k1,k2=k2,K=K,P=P,Q=Q,a=0,b=2,ij=FALSE,
                                                 dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                                 ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
            YijkWijl[j1,j2,k1,k2,i]<-E.YW.matrix(i=i,j1=j1,j2=j2,k1=k1,k2=k2,K=K,P=P,Q=Q,a=1,b=1,ij=FALSE,
                                                 dat=as.matrix(data),gam.d=gam.hat,beta.d=beta.hat,LRij.vec = LRij,
                                                 ibsMat = ibsMat,sigma.vector=sigma.vector,new.node=new.node,new.wt = new.wt,HIETAI = HIETAI)
          }}}}}
  
  
  T010<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=1,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T011<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=1,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T020<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=2,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T001<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=0,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T002<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=0,c=2,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T100<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T101<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T102<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=2,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T110<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=1,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T111<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=1,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T200<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)}) #E(eta_i^2)
  T201<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T202<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=2,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T210<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=1,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T300<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=3,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T301<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=3,b=0,c=1,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T400<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=4,b=0,c=0,d=0,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                       alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T0021<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=0,b=0,c=2,d=1,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                        alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T1021<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=1,b=0,c=2,d=1,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                        alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  T2021<-sapply(c(1:max.i),function(o){etai.abcd.vector(i=o,a=2,b=0,c=2,d=1,dat=as.matrix(data),P=P,Q=Q,new.wt=new.wt,new.node=node.mat[o,],
                                                        alpha.d=alpha.hat,c.d=c.hat,sigma.vector=sigma.vector,HIETAI=HIETAI)})
  
  d2QQ<-Edli2<-Sum.Edli.Edli<-matrix(0,nrow=K+P+Q+2,ncol=K+P+Q+2)
  
  d2Q1.c.sum<-0
  d2Q1.alpha.c.sum<-matrix(0,nrow=Q,ncol=1)
  d2Q1.gam.sum<-numeric(K)
  d2Q1.beta.sum<-matrix(0,nrow=P,ncol=P)
  d2Q1.gam.beta.sum<-matrix(0,nrow=K,ncol=P)
  d2Q1.alpha.sum<-matrix(0,nrow=Q,ncol=Q)
  d2Q1.alpha.c<-matrix(0,nrow=Q,ncol=1)
  for(i in c(1:max.i)){
    dat.ith<-as.matrix(data)[which(data[,1]==i),]
    if(is.matrix(dat.ith)==FALSE){
      Lij.vec<-dat.ith[3];Rij.vec<-dat.ith[4];
      DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
      Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
    }else{
      Lij.vec<-dat.ith[,3];Rij.vec<-dat.ith[,4]; 
      DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
      Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
    }
    dli.gam.k<-numeric(K);dli.beta<-numeric(P);dli.alpha<-numeric(Q)
    d2Q1.gam.i<-numeric(K)
    for(k in c(1:K)){d2Q1.gam.i[k]<-sum(sapply(c(1:length(Lij.vec)), function(o){gam.hat[k]^{-2}*(Yijk000[o,k,i]+Wijk000[o,k,i])}))}
    d2Q1.gam.sum<-d2Q1.gam.sum+d2Q1.gam.i
    
    d2Q1.beta.i<-matrix(0,nrow=P,ncol=P)
    for(p1 in c(1:P)){
      for(p2 in c(1:P)){
        d2Q1.beta.i[p1,p2]<-sum(sapply(c(1:length(Lij.vec)), function(o){
          Xij.vec[o,p1]*Xij.vec[o,p2]*exp(sum(Xij.vec[o,]*beta.hat))*T010[i]*
            ifelse(DeltaL.vec[o]==1,Lambda.Lij[i,o],Lambda.Rij[i,o])
        }))}}
    d2Q1.beta.sum<-d2Q1.beta.sum+d2Q1.beta.i
    
    d2Q1.gam.beta.i<-matrix(0,nrow=K,ncol=P)
    for(k in c(1:K)){
      for(p in c(1:P)){
        d2Q1.gam.beta.i[k,p]<-sum(sapply(c(1:length(Lij.vec)), function(o){
          Xij.vec[o,p]*exp(sum(Xij.vec[o,]*beta.hat))*T010[i]*
            ifelse(DeltaL.vec[o]==1,ibsMat[which(Lij.vec[o]==LRij),][k],ibsMat[which(Rij.vec[o]==LRij),][k])
        }))}}
    d2Q1.gam.beta.sum<-d2Q1.gam.beta.sum+d2Q1.gam.beta.i
    
    d2Q1.alpha.i<-matrix(0,nrow=Q,ncol=Q)
    for(q1 in c(1:Q)){for(q2 in c(1:Q)){d2Q1.alpha.i[q1,q2]<- max.m*Z.vec[q1]*Z.vec[q2]*T0021[i]}}
    d2Q1.alpha.sum<-d2Q1.alpha.sum+d2Q1.alpha.i
    
    d2Q1.alpha.c.i<-matrix(0,nrow=Q,ncol=1)
    for(q in c(1:Q)){d2Q1.alpha.c.i[q,1]<- max.m*Z.vec[q]*T1021[i]}
    d2Q1.alpha.c.sum<-d2Q1.alpha.c.sum+d2Q1.alpha.c.i
    
    d2Q1.c.sum<-d2Q1.c.sum+max.m*T2021[i]
    
    for(k in c(1:K)){dli.gam.k[k]<-sum(sapply(c(1:length(Lij.vec)), function(o){gam.hat[k]^{-1}*(Yijk000[o,k,i]+Wijk000[o,k,i])-
        exp(sum(Xij.vec[o,]*beta.hat))*T010[i]*
        ifelse(DeltaL.vec[o]==1,ibsMat[which(Lij.vec[o]==LRij),][k],ibsMat[which(Rij.vec[o]==LRij),][k])}))}
    
    for(p in c(1:P)){dli.beta[p]<-sum(sapply(c(1:length(Lij.vec)), function(o){Xij.vec[o,p]*(Yij000[i,o]+Wij000[i,o])-
        Xij.vec[o,p]*exp(sum(Xij.vec[o,]*beta.hat))*T010[i]*ifelse(DeltaL.vec[o]==1,Lambda.Lij[i,o],Lambda.Rij[i,o])}))}
    
    for(q in c(1:Q)){dli.alpha[q]<-Z.vec[q]*(length(Lij.vec)-max.m*T001[i])}
    dli.c<-length(Lij.vec)*T100[i]-max.m*T101[i]
    dli.sigma2<- -(2*sigma.hat^2)^(-1)+(2*sigma.hat^4)^(-1)*T200[i]
    
    dli<-as.matrix(c(dli.gam.k,dli.beta,dli.alpha,dli.c,dli.sigma2))
    Sum.Edli.Edli<-Sum.Edli.Edli+dli%*%t(dli)
  }
  
  d2Q2<-matrix(0,nrow=Q+1,ncol=Q+1)
  d2Q2[c(1:Q),c(1:Q)]<-d2Q1.alpha.sum; d2Q2[(Q+1),(Q+1)]<-d2Q1.c.sum
  d2Q2[c(1:Q),(Q+1)]<-d2Q1.alpha.c.sum;d2Q2[(Q+1),c(1:Q)]<-t(d2Q1.alpha.c.sum)
  d2Q3<- -max.i*(2*sigma.hat^4)^(-1)+sigma.hat^(-6)*sum(T200)
  
  diag(d2QQ[c(1:K),c(1:K)])<-d2Q1.gam.sum; d2QQ[c((K+1):(K+P)),c((K+1):(K+P))]<-d2Q1.beta.sum;
  d2QQ[c(1:K),c((K+1):(K+P))]<-d2Q1.gam.beta.sum; d2QQ[c((K+1):(K+P)),c(1:K)]<-t(d2Q1.gam.beta.sum)
  d2QQ[c((K+P+1):(K+P+Q+1)),c((K+P+1):(K+P+Q+1))]<-d2Q2; d2QQ[K+P+Q+2,K+P+Q+2]<-d2Q3
  
  G1<-matrix(0,nrow=K,ncol=K);G2<-matrix(0,nrow=P,ncol=P);G3<-matrix(0,nrow=Q,ncol=Q)
  G4<-matrix(0,nrow=P,ncol=K);G5<-matrix(0,nrow=P,ncol=Q);G6<-G7<-matrix(0,nrow=P,ncol=1);
  G8<-matrix(0,nrow=Q,ncol=K);G9<-matrix(0,nrow=Q,ncol=1);G10<-matrix(0,nrow=Q,ncol=1);G11<-G12<-matrix(0,nrow=1,ncol=K)
  
  dli2.c<-dli2.sigma<-dli2.c.sigma<-0
  for(i in c(1:max.i)){
    dat.ith<-as.matrix(data)[which(data[,1]==i),]
    if(is.matrix(dat.ith)==FALSE){
      Lij.vec<-dat.ith[3];Rij.vec<-dat.ith[4];
      DeltaL.vec<-dat.ith[5];DeltaI.vec<-dat.ith[6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[7])}else{Xij.vec<-t(as.matrix(dat.ith[(7:(7+P-1))]))}
      Z.vec  <-c(1,dat.ith[((7+P):(7+P+Q-2))])
    }else{
      Lij.vec<-dat.ith[,3];Rij.vec<-dat.ith[,4]; 
      DeltaL.vec<-dat.ith[,5];DeltaI.vec<-dat.ith[,6];
      if(P==1){Xij.vec<-as.matrix(dat.ith[,7])}else{Xij.vec<-dat.ith[,(7:(7+P-1))]}
      Z.vec  <-c(1,dat.ith[1,((7+P):(7+P+Q-2))])
    }
    
    dli2.gam.k<-matrix(0,nrow=K,ncol=K)
    for(u in c(1:K)){
      for(v in c(1:K)){
        dli2.gam.k[u,v]<-sum(sapply(c(1:length(Lij.vec)),function(o1){sapply(c(1:length(Lij.vec)),function(o2){
          if(o1!=o2){YY<-YijkYijl[o1,o2,u,v,i]
          }else if((o1==o2)&(u==v)){YY<- Yijk2[o1,u,i]
          }else if((o1==o2)&(u!=v)){YY<- -rijk[o1,u,i]*rijk[o1,v,i]*(Yij000[i,o1]-Yij2[i,o1])}
          if(o1!=o2){WW<-WijkWijl[o1,o2,u,v,i]
          }else if((o1==o2)&(u==v)){WW<- Wijk2[o1,u,i]
          }else if((o1==o2)&(u!=v)){WW<- -sijk[o1,u,i]*sijk[o1,v,i]*(Wij000[i,o1]-Wij2[i,o1])}
          
          (gam.hat[u]*gam.hat[v])^(-1)*(YY+YijkWijl[o1,o2,u,v,i]+YijkWijl[o2,o1,v,u,i]+WW)+exp(sum((Xij.vec[o1,]+Xij.vec[o2,])*beta.hat))*T020[i]*
            (ifelse(DeltaL.vec[o1]*DeltaL.vec[o2]==1,ibsMat[which(Lij.vec[o1]==LRij),][u]*ibsMat[which(Lij.vec[o2]==LRij),][v],0)+
               ifelse(DeltaL.vec[o1]*(1-DeltaL.vec[o2])==1,ibsMat[which(Lij.vec[o1]==LRij),][u]*ibsMat[which(Rij.vec[o2]==LRij),][v],0)+
               ifelse((1-DeltaL.vec[o1])*DeltaL.vec[o2]==1,ibsMat[which(Rij.vec[o1]==LRij),][u]*ibsMat[which(Lij.vec[o2]==LRij),][v],0)+
               ifelse((1-DeltaL.vec[o1])*(1-DeltaL.vec[o2])==1,ibsMat[which(Rij.vec[o1]==LRij),][u]*ibsMat[which(Rij.vec[o2]==LRij),][v],0))-
            gam.hat[u]^(-1)*(Yijk010[o1,u,i]+Wijk010[o1,u,i])*exp(sum(Xij.vec[o2,]*beta.hat))*
            (ifelse(DeltaL.vec[o2]==1,ibsMat[which(Lij.vec[o2]==LRij),][v],ibsMat[which(Rij.vec[o2]==LRij),][v]))-
            gam.hat[v]^(-1)*(Yijk010[o2,v,i]+Wijk010[o2,v,i])*exp(sum(Xij.vec[o1,]*beta.hat))*
            (ifelse(DeltaL.vec[o1]==1,ibsMat[which(Lij.vec[o1]==LRij),][u],ibsMat[which(Rij.vec[o1]==LRij),][u]))
        })}))
      }}
    G1<-G1+dli2.gam.k
    
    dli2.beta<-matrix(0,nrow=P,ncol=P)
    for(u in c(1:P)){
      for(v in c(1:P)){
        dli2.beta[u,v]<-sum(sapply(c(1:length(Lij.vec)),function(o1){sapply(c(1:length(Lij.vec)),function(o2){
          if(o1!=o2){YY<-YikYil[o1,o2,i]}else if(o1==o2){YY<-Yij2[i,o1]}
          if(o1!=o2){WW<-WikWil[o1,o2,i]}else if(o1==o2){WW<-Wij2[i,o1]}
          
          Xij.vec[o1,u]*Xij.vec[o2,v]*(YY+YikWil[o1,o2,i]+YikWil[o2,o1,i]+WW+
                                         exp(sum((Xij.vec[o1,]+Xij.vec[o2,])*beta.hat))*T020[i]*
                                         (ifelse(DeltaL.vec[o1]*DeltaL.vec[o2]==1,Lambda.Lij[i,o1]*Lambda.Lij[i,o2],0)+
                                            ifelse(DeltaL.vec[o1]*(1-DeltaL.vec[o2])==1,Lambda.Lij[i,o1]*Lambda.Rij[i,o2],0)+
                                            ifelse((1-DeltaL.vec[o1])*DeltaL.vec[o2]==1,Lambda.Rij[i,o1]*Lambda.Lij[i,o2],0)+
                                            ifelse((1-DeltaL.vec[o1])*(1-DeltaL.vec[o2])==1,Lambda.Rij[i,o1]*Lambda.Rij[i,o2],0))-
                                         (Yij010[i,o1]+Wij010[i,o1])*exp(sum(Xij.vec[o2,]*beta.hat))*
                                         (ifelse(DeltaL.vec[o2]==1,Lambda.Lij[i,o2],Lambda.Rij[i,o2]))-
                                         (Yij010[i,o2]+Wij010[i,o2])*exp(sum(Xij.vec[o1,]*beta.hat))*
                                         (ifelse(DeltaL.vec[o1]==1,Lambda.Lij[i,o1],Lambda.Rij[i,o1])))
        })}))
      }}
    G2<-G2+dli2.beta
    
    dli2.alpha<-matrix(0,nrow=Q,ncol=Q)
    for(u in c(1:Q)){
      for(v in c(1:Q)){
        dli2.alpha[u,v]<-(Z.vec[u]*Z.vec[v])*(length(Lij.vec)^2+max.m^2*T002[i]-2*length(Lij.vec)*max.m*T001[i])
      }}
    G3<-G3+dli2.alpha
    
    dli2.c    <-dli2.c + length(Lij.vec)^2*T200[i]+max.m^2*T202[i]-2*length(Lij.vec)*max.m*T201[i]
    dli2.sigma<-dli2.sigma + (4*sigma.hat^4)^(-1)+(4*sigma.hat^8)^(-1)*T400[i]-(2*sigma.hat^6)^(-1)*T200[i]
    
    dli2.beta.gammak<-matrix(0,nrow=P,ncol=K)
    for(u in c(1:P)){
      for(v in c(1:K)){
        dli2.beta.gammak[u,v]<-sum(sapply(c(1:length(Lij.vec)),function(o1){sapply(c(1:length(Lij.vec)),function(o2){
          if(o1!=o2){YY<-YikYijl[o1,o2,v,i]
          }else if(o1==o2){YY<- rijk[o1,v,i]*Yij2[i,o1]}
          if(o1!=o2){WW<-WikWijl[o1,o2,v,i]
          }else if(o1==o2){WW<- sijk[o1,v,i]*Wij2[i,o1]}
          
          Xij.vec[o1,u]*(gam.hat[v]^(-1)*(YY+YikWijl[o1,o2,v,i]+WikYijl[o1,o2,v,i]+WW)+exp(sum((Xij.vec[o1,]+Xij.vec[o2,])*beta.hat))*T020[i]*
                           (ifelse(DeltaL.vec[o1]*DeltaL.vec[o2]==1,Lambda.Lij[i,o1]*ibsMat[which(Lij.vec[o2]==LRij),][v],0)+
                              ifelse(DeltaL.vec[o1]*(1-DeltaL.vec[o2])==1,Lambda.Lij[i,o1]*ibsMat[which(Rij.vec[o2]==LRij),][v],0)+
                              ifelse((1-DeltaL.vec[o1])*DeltaL.vec[o2]==1,Lambda.Rij[i,o1]*ibsMat[which(Lij.vec[o2]==LRij),][v],0)+
                              ifelse((1-DeltaL.vec[o1])*(1-DeltaL.vec[o2])==1,Lambda.Rij[i,o1]*ibsMat[which(Rij.vec[o2]==LRij),][v],0))-
                           (Yij010[i,o1]+Wij010[i,o1])*exp(sum(Xij.vec[o2,]*beta.hat))*(ifelse(DeltaL.vec[o2]==1,
                                                                                               ibsMat[which(Lij.vec[o2]==LRij),][v],ibsMat[which(Rij.vec[o2]==LRij),][v]))-
                           gam.hat[v]^(-1)*(Yijk010[o2,v,i]+Wijk010[o2,v,i])*exp(sum(Xij.vec[o1,]*beta.hat))*
                           ifelse(DeltaL.vec[o1]==1,Lambda.Lij[i,o1],Lambda.Rij[i,o1])
          )})}))
      }}
    G4<-G4+dli2.beta.gammak
    
    dli2.beta.alpha<-matrix(0,nrow=P,ncol=Q)
    for(u in c(1:P)){
      for(v in c(1:Q)){
        dli2.beta.alpha[u,v]<-sum(sapply(c(1:length(Lij.vec)),function(o1){
          Xij.vec[o1,u]*Z.vec[v]*(length(Lij.vec)*(Yij000[i,o1]+Wij000[i,o1])+exp(sum(Xij.vec[o1,]*beta.hat))*(max.m*T011[i]-length(Lij.vec)*T010[i])*
                                    ifelse(DeltaL.vec[o1]==1,Lambda.Lij[i,o1],Lambda.Rij[i,o1])-max.m*(Yij001[i,o1]+Wij001[i,o1]))
        }))}}
    G5<-G5+dli2.beta.alpha
    
    dli2.beta.c<-matrix(0,nrow=P,ncol=1)
    for(u in c(1:P)){
      dli2.beta.c[u,1]<-sum(sapply(c(1:length(Lij.vec)),function(o1){
        Xij.vec[o1,u]*(length(Lij.vec)*(Yij100[i,o1]+Wij100[i,o1])+exp(sum(Xij.vec[o1,]*beta.hat))*(max.m*T111[i]-length(Lij.vec)*T110[i])*
                         ifelse(DeltaL.vec[o1]==1,Lambda.Lij[i,o1],Lambda.Rij[i,o1])-max.m*(Yij101[i,o1]+Wij101[i,o1]))
      }))}
    G6<-G6+dli2.beta.c
    
    dli2.beta.sigma<-matrix(0,nrow=P,ncol=1)
    for(u in c(1:P)){
      dli2.beta.sigma[u,1]<-sum(sapply(c(1:length(Lij.vec)),function(o1){
        Xij.vec[o1,u]*(-(2*sigma.hat^2)^(-1)*(Yij000[i,o1]+Wij000[i,o1])+(2*sigma.hat^2)^(-1)*exp(sum(Xij.vec[o1,]*beta.hat))*(
          T010[i]-(sigma.hat)^(-2)*T210[i])*ifelse(DeltaL.vec[o1]==1,Lambda.Lij[i,o1],Lambda.Rij[i,o1])+
            (2*sigma.hat^4)^(-1)*(Yij200[i,o1]+Wij200[i,o1]))
      }))}
    G7<-G7+dli2.beta.sigma
    
    dli2.alpha.gam.k<-matrix(0,nrow=Q,ncol=K)
    for(u in c(1:Q)){
      for(v in c(1:K)){
        dli2.alpha.gam.k[u,v]<-sum(sapply(c(1:length(Lij.vec)),function(o1){
          Z.vec[u]*(length(Lij.vec)*gam.hat[v]^(-1)*(Yijk000[o1,v,i]+Wijk000[o1,v,i])+exp(sum(Xij.vec[o1,]*beta.hat))*
                      (max.m*T011[i]-length(Lij.vec)*T010[i])*
                      ifelse(DeltaL.vec[o1]==1,ibsMat[which(Lij.vec[o1]==LRij),][v],ibsMat[which(Rij.vec[o1]==LRij),][v])-
                      max.m*gam.hat[v]^(-1)*(Yijk001[o1,v,i]+Wijk001[o1,v,i]))
        }))}}
    G8<-G8+dli2.alpha.gam.k
    
    dli2.alpha.c<-matrix(0,nrow=Q,ncol=1)
    for(u in c(1:Q)){dli2.alpha.c[u,1]<-Z.vec[u]*(length(Lij.vec)^2*T100[i]+max.m^2*T102[i]-2*max.m*length(Lij.vec)*T101[i])}
    G9<-G9+dli2.alpha.c
    
    dli2.alpha.sigma<-matrix(0,nrow=Q,ncol=1)
    for(u in c(1:Q)){dli2.alpha.sigma[u,1]<-Z.vec[u]*(-(2*sigma.hat^2)^(-1)*(length(Lij.vec)-max.m*T001[i])+(2*sigma.hat^4)^(-1)*
                                                        (length(Lij.vec)*T200[i]-max.m*T201[i]))}
    G10<-G10+dli2.alpha.sigma
    
    dli2.c.gam.k<-matrix(0,nrow=1,ncol=K)
    for(v in c(1:K)){dli2.c.gam.k[1,v]<-sum(sapply(c(1:length(Lij.vec)),function(o2){
      length(Lij.vec)*gam.hat[v]^(-1)*(Yijk100[o2,v,i]+Wijk100[o2,v,i])+exp(sum(Xij.vec[o2,]*beta.hat))*(
        max.m*T111[i]-length(Lij.vec)*T110[i])*ifelse(DeltaL.vec[o2]==1,ibsMat[which(Lij.vec[o2]==LRij),][v],ibsMat[which(Rij.vec[o2]==LRij),][v])-
        max.m*gam.hat[v]^(-1)*(Yijk101[o2,v,i]+Wijk101[o2,v,i])
    }))}
    G11<-G11+dli2.c.gam.k
    
    dli2.c.sigma<-dli2.c.sigma-(2*sigma.hat^2)^(-1)*(length(Lij.vec)*T100[i]-max.m*T101[i])+(2*sigma.hat^4)^(-1)*(length(Lij.vec)*T300[i]-max.m*T301[i])
    
    dli2.sigma.gam.k<-matrix(0,nrow=1,ncol=K)
    for(v in c(1:K)){dli2.sigma.gam.k[1,v]<-sum(sapply(c(1:length(Lij.vec)),function(o2){
      -(gam.hat[v]*(2*sigma.hat^2))^(-1)*(Yijk000[o2,v,i]+Wijk000[o2,v,i])+(2*sigma.hat^2)^(-1)*exp(sum(Xij.vec[o2,]*beta.hat))*
        (T010[i]-sigma.hat^(-2)*T210[i])*ifelse(DeltaL.vec[o2]==1,ibsMat[which(Lij.vec[o2]==LRij),][v],ibsMat[which(Rij.vec[o2]==LRij),][v])+
        (gam.hat[v]*2*sigma.hat^4)^(-1)*(Yijk200[o2,v,i]+Wijk200[o2,v,i]) 
    }))}
    G12<-G12+dli2.sigma.gam.k
  }
  
  Edli2[c(1:K),c(1:K)]<-G1;Edli2[c((K+1):(K+P)),c((K+1):(K+P))]<-G2;Edli2[c((K+P+1):(K+P+Q)),c((K+P+1):(K+P+Q))]<-G3;
  Edli2[(K+P+Q+1),(K+P+Q+1)]<-dli2.c;Edli2[(K+P+Q+2),(K+P+Q+2)]<-dli2.sigma;Edli2[c(1:K),c((K+1):(K+P))]<-t(G4);
  Edli2[c((K+1):(K+P)),c((K+P+1):(K+P+Q))]<-G5;Edli2[c((K+1):(K+P)),(K+P+Q+1)]<-G6;Edli2[c((K+1):(K+P)),(K+P+Q+2)]<-G7;
  Edli2[c(1:K),c((K+P+1):(K+P+Q))]<-t(G8);Edli2[c((K+P+1):(K+P+Q)),(K+P+Q+1)]<-G9;Edli2[c((K+P+1):(K+P+Q)),(K+P+Q+2)]<-G10;
  Edli2[c(1:K),K+P+Q+1]<-t(G11);Edli2[(K+P+Q+1),(K+P+Q+2)]<-dli2.c.sigma;Edli2[c(1:K),K+P+Q+2]<-t(G12)
  
  for(i in c(1:(K+P+Q+2))){
    for(j in c(1:(K+P+Q+2))){
      if(i>j){Edli2[i,j]<-Edli2[j,i]}}}
  
  sd.array<-sqrt(diag(solve(d2QQ-Edli2+Sum.Edli.Edli)))
  beta.se<-sd.array[c((K+1):(K+P))]
  alpha.se<-sd.array[c((K+P+1):(K+P+Q))]
  c.se<-sd.array[(K+P+Q+1)]
  sigma.se<-sd.array[(K+P+Q+2)]*(2*sigma.hat)^(-1)
  
  list(loglik= logLike, gam.hat=gam.hat,
       alpha.hat=alpha.hat, kappa.hat=c.hat, beta.hat=beta.hat, sigma.hat=sigma.hat,
       alpha.hat.se=alpha.se, kappa.hat.se=c.se, beta.hat.se=beta.se, sigma.hat.se=sigma.se)
}


