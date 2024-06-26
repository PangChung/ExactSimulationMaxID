### adaptive rejection sampler ###
# logarithmic target function #

library(ars)
logends<-function(x,y,parR){
  val <- c()
  alpha = parR[1];beta = parR[2];nu=parR[3]
  ind <- !is.na(x)
  val[!ind] <- NA
  r<-x;r[ind]<-exp(x[ind])
  y = y/sqrt(nu)
  if(beta==0){
    val[ind] <- -1/2*(y/r[ind])^2 - (alpha+1)*log(r[ind])+log(alpha)
  }else{
    val[ind] <- -1/2*(y/r[ind])^2-alpha/beta*(r[ind]^beta-1)+log(beta*r[ind]^(-beta-1)+alpha*r[ind]^(-1))
  }
  return(val)
}

## the derivative of the logarithmic target function ##
grad_logends<-function(x,y,parR){
  alpha <- parR[1]
  beta <- parR[2];nu=par[3]
  val <- c()
  val[is.na(x)] <- NA
  ind <- !is.na(x)
  r<-x;r[ind]<-exp(x[ind])
  y = y/sqrt(nu)
  if(beta==0){
    val[ind] <- y^2*r[ind]^(-2) - (alpha+1)
  }else{
    val[ind] <-  y^2*r[ind]^(-2) - alpha*r[ind]^(beta) - (beta*(beta+1)*r[ind]^(-beta-1)+alpha*r[ind]^(-1))/(beta*r[ind]^(-beta-1)+alpha*r[ind]^(-1))
  }
  return(val)
}

ars.y <- function(y,idx,parR,Sigma){
  init <- seq(-4,4,length=6)
  log.r <- ars(n=1,logends,grad_logends,x=init,y=y,m=length(init),parR=parR)
  Y=rcondmv(y=y,r=exp(log.r),idx=idx,Sigma=Sigma,parR=parR)
  return(Y)
}

#ars(1,logends,grad_logends,x=c(-3,-1,0,1,3,4),y=0.5,parR=c(1,1))

#### some extra functions to use for our specific model ####
rcondmv <- function(y,r,idx=1,Sigma,parR,ncores=1){
  D <- nrow(Sigma)
  Z <- rep(NA,D)
  Z[idx] <- y
   A <- Sigma[-idx,-idx]-Sigma[-idx,idx]%*%t(Sigma[-idx,idx])/Sigma[idx,idx]
   A <- t(chol(A))
   Z[-idx] <- r*A %*% matrix(rnorm(D-1)) + Sigma[-idx,idx]/Sigma[idx,idx]*y
   return(Z)
}

### mcmc alogrithm ###
mcmc.y <- function(y,idx,parR,Sigma,sd=1,N=10^3){
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
  } 
  r = exp(log.r[count])
  ## simulate X1 | X2=y, R
  Y=rcondmv(y=y,r=r,idx=idx,Sigma=Sigma,parR=parR)
  return(Y)
}


## Algorithm to sample using MH algorithm or adaptive rejection sampling ##
mh <- function(n=1,ars=F,parR,Sigma,N=10^3){
  counter = 1
  z = rexp(1)
  y = qG(exp(-z),parR = parR)
  D = nrow(Sigma)
  if(!ars){
    Y = mcmc.y(y = y,idx=1,parR=parR,Sigma=Sigma,sd=0.5,N=N)
  }else{
    Y = ars.y(y = y,idx=1,parR=parR,Sigma=Sigma)
  }
  for(j in 2:D){
    z = rexp(1)
    y = qG(exp(-z),parR = parR)
    while(y > Y[j]){
      if(!ars){
        Y.new = mcmc.y(y = y,idx=j,parR=parR,Sigma=Sigma,sd=0.5,N=N)
      }else{
        Y.new = ars.y(y = y,idx=j,parR=parR,Sigma=Sigma)
      }
      if(!any(Y.new[1:(j-1)] >= Y[1:(j-1)])){
        Y = pmax(Y.new,Y)
      }
      e = rexp(1)
      z = z + e
      y = qG(exp(-z),parR=parR)
      counter = counter + 1
    }
  }
  return(Y)
}

# mh <- function(n=1,ars=F,parR,Sigma,N=10^3){
#     D = nrow(Sigma)
#     idx.loc = cbind(1:n,1)
#     idx.loc.sum = lapply(1:D,function(id){which(idx.loc[,2]==id)})
#     z = rexp(n)
#     y = qG(exp(-z),parR = parR)
#     func <- function(idx.j,j){
#         m.idx.j = length(idx.j)
#         if(m.idx.j > 0 & j<=n){
#             if(ars){
#                 val = lapply(idx.j,function(i){ars.y(y = y[i],idx=j,parR=parR,Sigma=Sigma)})
#             } else {
#                 val = lapply(idx.j,function(i){mcmc.y(y = y[i],idx=j,parR=parR,Sigma=Sigma,sd=0.5,N=N)})
#             }
#             return(do.call(rbind,val))
#         }
#         return(NULL)
#     }

#     func.compare <- function(idx,j){
#         val = !any(Y[idx,1:(j-1)] < Y.temp[idx,1:(j-1)]) 
#         if(val){
#             return(idx)
#         }else{
#             return(NULL)
#         }
#     }
    
#     func.pmax <- function(idx){
#         return(pmax(Y[idx,],Y.temp[idx,]))
#     }

#     Y = mapply(func,idx.j=idx.loc.sum,j=1:D,SIMPLIFY=FALSE)
#     Y = do.call(rbind,Y) 
#     idx.loc[,2] = 2
#     idx.loc.sum = lapply(1:D,function(id){which(idx.loc[,2]==id)})
#     idx.new <- idx.pass <-  rep(TRUE,n)
#     idx.finish <- idx.old <- rep(FALSE,n)
#     count = 1
#     while(any(!idx.finish)){
#         if(sum(idx.old)!=0){
#             z[idx.old] = z[idx.old ] + rexp(sum(idx.old))
#             y[idx.old] = qG(exp(-z[idx.old]),parR)
#         }
#         if(sum(idx.new) != 0){
#             z[idx.new] = rexp(sum(idx.new))
#             y[idx.new] = qG(exp(-z[idx.new]),parR)
#         }
#         idx.pass[!idx.finish] <- y[!idx.finish] > Y[idx.loc[!idx.finish,,drop=FALSE]]
#         idx.pass[idx.finish] = FALSE
#         idx.go.ahead = !idx.pass & !idx.finish
#         if(sum(idx.go.ahead) > 0){
#             idx.loc[idx.go.ahead,2] = idx.loc[idx.go.ahead,2] +  1
#             idx.new[idx.go.ahead] = TRUE
#         }
#         if(sum(idx.pass)!=0){
#             idx.new[idx.pass] = FALSE
#             idx.loc.sum.simu = lapply(1:n,function(id){which(idx.loc[idx.pass,2]==id)})
#             Y.temp.0 = mapply(func,idx.j=idx.loc.sum.simu,j=1:n,SIMPLIFY=FALSE)
#             Y.temp = Y
#             idx.temp = which(idx.pass)[unlist(idx.loc.sum.simu)]
#             Y.temp[idx.temp,] = do.call(rbind,Y.temp.0)
#             idx.temp = unlist(mapply(func.compare,idx=which(idx.pass),j=idx.loc[idx.pass,2],SIMPLIFY=FALSE))
#             if(length(idx.temp)!=0){
#                 Y[idx.temp,] = do.call(rbind,lapply(idx.temp,func.pmax))
#                 idx.loc[idx.temp,2] = idx.loc[idx.temp,2] + 1
#                 idx.new[idx.temp] = TRUE
#             }
#         }
#         idx.finish = idx.loc[,2] > D
#         idx.new = idx.new & !idx.finish
#         idx.old = !idx.new & !idx.finish
#         count = count + 1
#         # print(paste0(c(count,sum(idx.finish),n),collapse = "/"))
#     }
#     return(Y)
# }

emp.pair <- function(k,pair,z=1,XDAT,parR,cor.mat){# index of the pairs
  set.seed(19873436)
  d <- ncol(pair)
  dat0 <- XDAT[,pair[k,]]
  dat.max <- apply(dat0,1,max,na.rm=T)
  n <- sum(dat.max<z)
  pf<-ecdf(dat.max)
  p <- pf(z)
  v <- pmin(d,pmax(1,-z*log(p)))
  sd <- 2*sqrt((1-p)/n*z^2/p)
  theta <- Theta.dimD(z=z,parR = parR,cor.mat[pair[k,],pair[k,]])
  return(cbind(v,sd,theta))
}



