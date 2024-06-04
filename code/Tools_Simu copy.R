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

ars.y <- function(y,idx,pars,coord,reg,reg.t){
  init <- seq(-4,4,length=6)
  #log.r <- ars(n=n,logends,grad_logends,x=init,m=length(init),y=y,parR=pars$parR)
  #init <- sort(c(init,log.r))
  log.r <- ars(n=1,logends,grad_logends,x=init,y=y,m=length(init),parR=pars$parR)
  Y=rcondmv(y=y,r=exp(log.r),idx=idx,coord=coord,parR=pars$parR,parGauss=pars$parGauss,reg=reg,reg.t=reg.t)
  return(Y)
}

#ars(1,logends,grad_logends,x=c(-3,-1,0,1,3,4),y=0.5,parR=c(1,1))

#### some extra functions to use for our specific model ####
rcondmv <- function(y,r,idx=1,coord, parR, parGauss,reg=NULL,reg.t=NULL,ncores=1,type="GS"){
  D <- nrow(coord)
  Z <- rep(NA,D)
  Z[idx] <- y
  pairs <- combn(1:D,2)
  Sigma <- diag(1,D)
  if(is.null(parGauss$nu)){
    fun <- function(pair){
      if(parGauss$type %in% c(1,2,4)){
        h = coord[pair[1],pair[2]]
      }else{h = abs(coord[pair[1],]-coord[pair[2],])}
      reg.new=NULL;reg.t.new=NULL
      if(!is.null(reg)){reg.new=reg[pair,]}
      if(!is.null(reg.t)){reg.t.new=reg.t}
      rho <- rho.func(h=h,r=NULL,parGauss =parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
      return(rho)
    }
      Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
      Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
  } else{	
        fun <- function(pair){
          if(parGauss$type %in% c(1,2,4)){
            h = coord[pair[1],pair[2]]
          }else{h = abs(coord[pair[1],]-coord[pair[2],])}
          if(!is.null(reg)){reg.new=reg[pair,]}else{reg.new=NULL}
          if(!is.null(reg.t)){reg.t.new=reg.t}else{reg.t.new=NULL}
          rho <- rho.func(h=h,r=r,parGauss=parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
          return(rho)
        }
        Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
        Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
  }
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

mcmc.y <- function(y,idx,pars,coord,reg,reg.t,sd=1,N=10^3,type="GS"){
  if(type=="GS"){
  log.r = rnorm(1,mean=log(y),sd=sd) # proposal for R (hope this would be fine)
  logaccp = NULL
  ### algos to simulate the RW
  count = 1
  while(count <= N){
    ## simulate R|X2 = y
    log.r[count+1] = rnorm(1,mean = log.r[count],sd=sd)
    logaccp[count] <- min(0,logends(log.r[count+1],y,pars$parR) - logends(log.r[count],y,pars$parR) )
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
      logaccp[count] <- min(0,logends(log.r[count+1],y,pars$parR,type=type) - logends(log.r[count],y,pars$parR,type=type) )
      #message(count," :",log.r[count+1])
      if(runif(1) > exp(logaccp[count])){
        log.r[count+1] = log.r[count]
      }
      count  = count + 1
    } 
    r=log.r[count] }
  ## simulate X1 | X2=y, R
  Y=rcondmv(y=y,r=r,idx=idx,coord=coord,parR=pars$parR,parGauss=pars$parGauss,reg=reg,reg.t=reg.t,type=type)
  return(Y)
}


#### Algorithm to sample using MH algorithm or adaptive rejection sampling ####
mh <- function(n=1,ars=F,reg,reg.t,pars,coord,N=10^3,type="GS"){
  counter = 1
  z = rexp(1)
  y = qG(exp(-z),parR = pars$parR)
  D = nrow(coord)
  if(!ars){
    Y = mcmc.y(y = y,idx=1,pars=pars,coord=coord,reg=reg,reg.t=reg.t,sd=0.1,N=N,type = type)
  }else{
    Y = ars.y(y = y,idx=1,pars=pars,coord=coord,reg=reg,reg.t=reg.t)
  }
  for(j in 2:D){
    z = rexp(1)
    y = qG(exp(-z),parR = pars$parR,type = type)
    while(y > Y[j]){
      if(!ars){
        Y.new = mcmc.y(y = y,idx=j,pars=pars,coord=coord,reg=reg,reg.t=reg.t,sd=0.1,N=N,type=type)
      }else{
        Y.new = ars.y(y = y,idx=j,pars=pars,coord=coord,reg=reg,reg.t=reg.t)
      }
      if(!any(Y.new[1:(j-1)] >= Y[1:(j-1)])){
        Y = pmax(Y.new,Y)
      }
      e = rexp(1)
      z = z + e
      y = qG(exp(-z),parR = pars$parR,type = type)
      counter = counter + 1
    }
  }
  return(list(Y=Y,counter=counter))
}


## functions that uesd in the bounded method ##
upperG <- function (x, parR, d=49,log = FALSE){
  g <- function(xi) {
    if (xi <= 0) {
      return(Inf)
    }
    fun <- function(r, parR,di) {
      return(exp(dgamma(r^2,shape = di/2,scale = 2,log = T) + log(2*r) + dF(xi/r, parR, log = TRUE)))
    }
    val <- integrate(fun, lower = 0, upper = Inf, parR = parR, di=d,
                     rel.tol = 10^(-8), stop.on.error = FALSE)$value
    return(val)
  }
  val <- c()
  I <- !is.na(x)
  val[I] <- apply(matrix(x[I], ncol = 1), 1, g)
  val[!I] <- NA
  logval <- log(val)
  if (log) {
    return(logval)
  }
  else {
    return(exp(logval))
  }
}

upperGinv <- function (y, parR,d=49,log = FALSE) 
{
  alpha <- parR[1]
  beta <- parR[2]
  logval <- c()
  for (i in 1:length(y)) {
    if (!is.na(y[i]) & y[i] > 0) {
      fun <- function(x) {
        return(log(y[i]) - upperG(x = exp(x), parR = parR,d=d,
                                  log = TRUE))
      }
      logval[i] <- uniroot(f = fun, interval = c(-3, 3), 
                           extendInt = "yes")$root
    }
  }
  if (log) {
    return(logval)
  }
  else {
    return(exp(logval))
  }
}
### compute the empirical extremal coefficients ###
#D <- nrow(coord)
#z <- qfrechet(c(0.05,0.25,0.5,0.75,0.95))
#pairs2 <- t(combn(49,2))
#need coord,XDAT,pars

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
### effective truncated sampler algorithm ###

rG <- function(parGauss,parR,coord,reg,reg.t,ncores=1){
  D <- nrow(coord)
  Z <- rep(0,D)
  pairs <- combn(1:D,2)
  sigma.func <- function(r){
  Sigma <- diag(1,D)
  if(is.null(parGauss$nu)){
    fun <- function(pair){
      if(parGauss$type %in% c(1,2,4)){
        h = coord[pair[1],pair[2]]
      }else{h = abs(coord[pair[1],]-coord[pair[2],])}
      reg.new=NULL;reg.t.new=NULL
      if(!is.null(reg)){reg.new=reg[pair,]}
      if(!is.null(reg.t)){reg.t.new=reg.t}
      rho <- rho.func(h=h,r=NULL,parGauss =parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
      return(rho)
    }
    if(is.null(reg.t)){
      Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
      Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
    }else{
      Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
      Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
    }
  } else{	
    fun <- function(pair){
      if(parGauss$type %in% c(1,2,4)){
        h = coord[pair[1],pair[2]]
      }else{h = abs(coord[pair[1],]-coord[pair[2],])}
      if(!is.null(reg)){reg.new=reg[pair,]}else{reg.new=NULL}
      if(!is.null(reg.t)){reg.t.new=reg.t}else{reg.t.new=NULL}
      rho <- rho.func(h=h,r=r,parGauss=parGauss,reg=reg.new,reg.t=reg.t.new,mat=F)
      return(rho)
    }
    Sigma[t(pairs)]<-mcmapply(fun, pair=as.list(as.data.frame(pairs)),mc.cores = ncores)
    Sigma[t(pairs)[,2:1]]<-Sigma[t(pairs)]
  }
   return(Sigma)
  }
  r <- rexp(1); R <- upperGinv(exp(-r),parR=parR,d = D)
  A <- t(chol(sigma.func(R)))
  A.max <- max(apply(A,1,sum))
  S <- rnorm(D) ; S <- S/sqrt(sum(S^2))
  Z = R*A%*%S
  while(any(R*A.max > Z)){
    r <- r + rexp(1)
    print(paste(R*A.max,min(Z)))
    R <- upperGinv(r,parR=parR,d = D)
    A <- t(chol(sigma.func(R)))
    A.max <- max(apply(A,1,sum))
    S <- rnorm(D) ; S <- S/sqrt(sum(S^2))
    Z = pmax(R*A%*%S,Z)
  }
  return(Z)
}






