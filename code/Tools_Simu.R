### adaptive rejection sampler ###
# logarithmic target function #

library(ars)
logends<-function(x,y,parR,type="GS"){
  val <- c()
  alpha = parR[1];beta = parR[2]
  ind <- !is.na(x)
  val[!ind] <- NA
  r<-x;r[ind]<-exp(x[ind])
  if(type=="GS"){
  if(beta==0){
    val[ind] <- -1/2*(y/r[ind])^2 - (alpha+1)*log(r[ind])+log(alpha)
  }else{
    val[ind] <- -1/2*(y/r[ind])^2-alpha/beta*(r[ind]^beta-1)+log(beta*r[ind]^(-beta-1)+alpha*r[ind]^(-1))
  }}else{
    val[ind] <- -1/2*(y-x[ind])^2 + dF_GL(x[ind],parR,log=TRUE)
  }
  return(val)
}

## the derivative of the logarithmic target function ##
grad_logends<-function(x,y,parR){
  alpha <- parR[1]
  beta <- parR[2]
  val <- c()
  val[is.na(x)] <- NA
  ind <- !is.na(x)
  r<-x;r[ind]<-exp(x[ind])
  if(beta==0){
    val[ind] <- y^2*r[ind]^(-2) - (alpha+1)
  }else{
    val[ind] <-  y^2*r[ind]^(-2) - alpha*r[ind]^(beta) - (beta*(beta+1)*r[ind]^(-beta-1)+alpha*r[ind]^(-1))/(beta*r[ind]^(-beta-1)+alpha*r[ind]^(-1))
  }
  return(val)
}

ars.y <- function(y,idx,parR,coord){
  init <- seq(-4,4,length=6)
  log.r <- ars(n=1,logends,grad_logends,x=init,y=y,m=length(init),parR=parR)
  Y=rcondmv(y=y,r=exp(log.r),idx=idx,coord=coord,parR=parR)
  return(Y)
}

#ars(1,logends,grad_logends,x=c(-3,-1,0,1,3,4),y=0.5,parR=c(1,1))

#### some extra functions to use for our specific model ####
rcondmv <- function(y,r,idx=1,coord,parR,ncores=1,type="GS"){
  D <- nrow(coord)
  Z <- rep(NA,D)
  Z[idx] <- y
   Sigma=coord
   A <- Sigma[-idx,-idx]-Sigma[-idx,idx]%*%t(Sigma[-idx,idx])/Sigma[idx,idx]
   A <- t(chol(A))
   if(type=="GS"){
   Z[-idx] <- r*A %*% matrix(rnorm(D-1)) + Sigma[-idx,idx]/Sigma[idx,idx]*y
   }else{
     Z[-idx] <- A %*% matrix(rnorm(D-1)) + Sigma[-idx,idx]/Sigma[idx,idx]*(y-r)+r
   }
   return(Z)
}

### mcmc alogrithm ###

mcmc.y <- function(y,idx,parR,coord,sd=1,N=10^3,type="GS"){
  if(type=="GS"){
  log.r = rnorm(1,mean=log(y),sd=sd) # proposal for R (hope this would be fine)
  logaccp = NULL
  ### algos to simulate the RW
  count = 1
  while(count <= N){
    ## simulate R|X2 = y
    log.r[count+1] = rnorm(1,mean = log.r[count],sd=sd)
    logaccp[count] <- min(0,logends(log.r[count+1],y,parR) - logends(log.r[count],y,parR) )
    #message(count," :",log.r[count+1])
    if(runif(1) > exp(logaccp[count])){
      log.r[count+1] = log.r[count]
    }
    count  = count + 1
    #   ## plot r per 1000 iterations   
    #   if(count %% 1000 == 0){
    #     par(mfrow=c(1,2),mar=c(4,4,4,1),cex=0.5,cex.main=1,cex.lab=1,mgp=c(3,1,0))
    #     plot(1:count,exp(log.r),xlab="Iterations",ylab="log R",main="",pch=16,type="l")
    #     plot(1:(count-1),exp(logaccp),xlab="Iterations",ylab="log accp",main="",pch=16,type="l")
    #   }
  } 
  r = exp(log.r[count]) }else{
    log.r = rnorm(1,mean=y,sd=sd) # proposal for R (hope this would be fine)
    logaccp = NULL
    ### algos to simulate the RW
    count = 1
    while(count <= N){
      ## simulate R|X2 = y
      log.r[count+1] = rnorm(1,mean = log.r[count],sd=sd)
      logaccp[count] <- min(0,logends(log.r[count+1],y,parR,type=type) - logends(log.r[count],y,parR,type=type) )
      #message(count," :",log.r[count+1])
      if(runif(1) > exp(logaccp[count])){
        log.r[count+1] = log.r[count]
      }
      count  = count + 1
    } 
    r=log.r[count] }
  ## simulate X1 | X2=y, R
  Y=rcondmv(y=y,r=r,idx=idx,coord=coord,parR=parR,type=type)
  return(Y)
}


#### Algorithm to sample using MH algorithm or adaptive rejection sampling ####
mh <- function(n=1,ars=F,parR,coord,N=10^3,type="GS"){
  counter = 1
  z = rexp(1)
  y = qG(exp(-z),parR = parR,type=type)
  D = nrow(coord)
  if(!ars){
    Y = mcmc.y(y = y,idx=1,parR=parR,coord=coord,sd=0.01,N=N,type = type)
  }else{
    Y = ars.y(y = y,idx=1,parR=parR,coord=coord)
  }
  for(j in 2:D){
    z = rexp(1)
    y = qG(exp(-z),parR = parR,type = type)
    while(y > Y[j]){
      if(!ars){
        Y.new = mcmc.y(y = y,idx=j,parR=parR,coord=coord,sd=0.01,N=N,type=type)
      }else{
        Y.new = ars.y(y = y,idx=j,parR=parR,coord=coord)
      }
      if(!any(Y.new[1:(j-1)] >= Y[1:(j-1)])){
        Y = pmax(Y.new,Y)
      }
      e = rexp(1)
      z = z + e
      y = qG(exp(-z),parR =parR,type = type)
      counter = counter + 1
    }
  }
  return(list(Y=Y,counter=counter))
}


emp.pair <- function(k,pair,XDAT,pars,reg,reg.t){# index of the pairs
  set.seed(19873436)
  d <- ncol(pair)
  dat0 <- XDAT[,pair[k,]]
  A0 <- apply(dat0,1,function(x){sum(!is.na(x))==length(x)}) # time without missing value for the pair
  n0 <- sum(A0)
  if (n0<30) return(matrix(NA,length(z),ncol=3))
  dat.max <- apply(dat0[A0,],1,max,na.rm=T)
  n <- sapply(z,function(x){sum(dat.max<x)},simplify = T)
  pf<-ecdf(dat.max)
  p <- pf(z)
  A1 <- (p == 0)
  p[A1]<-NA
  v = rep(NA,length(z))
  v[!A1] <- pmin(d,pmax(1,-z[!A1]*log(p[!A1])))
  sd <- 2*sqrt((1-p)/n*z^2/p)
  if(d==2){h=coord[pair[k,1],pair[k,2]]}else{h=coord[pair[k,],pair[k,]]}
  theta <- Theta.dimD(z=z,h=h,parR = pars$parR,parGauss = pars$parGauss,reg = reg[pair[k,],],reg.t=reg.t)
  return(cbind(v,sd,theta))
}


message("Input the covariance matrix as coord")


