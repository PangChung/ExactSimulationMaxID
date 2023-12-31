args <- commandArgs(trailingOnly=TRUE)
for(arg in args) {
  eval(parse(text = arg))
}

setwd("~/project/exact_simu_max_id_project/")
source("Tools_MaxID.R");source("Tools_Simu.R")
ncores=40
library(parallel)
library(mvtnorm)
library(evd)

d = 7
D = d^2 # number of sites in the space
N = 100   # number of iterations when sample Y 
n = 50
LAMBDA.T = TRUE # whether the exponent measure is infinite or not on the whole space
parR = c(1,1) #alpha, beta
parGauss = c(1,1,0,1) #nu, lambdas
pars <- get.par(c(log(parR),parGauss),4)
x.coord <-  y.coord <- c(1:d)/(d+1) #grids 
cutoff <- 3/(7+1) # six order neighbors
coord = as.matrix(expand.grid(x.coord,y.coord))
coord = as.matrix(dist(coord))
reg=cbind(1,2*pnorm(coord[,1],0.5,0.25)-1)
reg.t = 0

u.approx <- apply(mcmapply(rmaxidspat,n=rep(1,n),MoreArgs = list(coord=coord,parR=pars$parR,parGauss = pars$parGauss,reg=reg,reg.t=reg.t),SIMPLIFY = T,mc.cores = ncores),1,pG,parR=pars$parR)

u.mh <- apply(mcmapply(mh,n=1:n,SIMPLIFY = T,mc.cores = ncores,MoreArgs = list(reg=reg,reg.t=reg.t,pars=pars)),1,pG,parR=pars$parR)

u.ars <- apply(mcmapply(mh,n=1:n,MoreArgs = list(ars=T,reg=reg,reg.t=reg.t,pars=pars),SIMPLIFY = T,mc.cores = ncores),1,pG,parR=pars$parR)


fixed <- c(F,F,F,F,T,F) ## corresponding to c(alpha, beta, lambda.0,lambda.1,lambda.2,nu)
init <- c(pars$parR,pars$parGauss$lambda,pars$parGauss$lambda.t,pars$parGauss$nu)
init[!fixed] = init[!fixed] + runif(sum(!fixed),0.1,0.3)
init[1:2] <- log(init[1:2])

### Fit the model. It will take hours. ###
message("Fit the model for ARS")
fit.ars <- fit.pw.parallel(init=init,datU = u.ars,coord=coord,reg=reg,reg.t=NULL,cutoff=cutoff,proppairs= 1,fixed=fixed,optim =T, hessian=F,sandwich=F,eps = 10^(-2), print.par.file=NULL,ncores=ncores,fit.load=F, fit.save=F,fit.file=NULL)

message("Fit the model for Metropolis Hastings")
fit.mh <- fit.ars <- fit.pw.parallel(init=init,datU = u.mh,coord=coord,reg=reg,reg.t=NULL,cutoff=cutoff,proppairs= 1,fixed=fixed,optim =T, hessian=F,sandwich=F,eps = 10^(-2), print.par.file=NULL,ncores=ncores,fit.load=F, fit.save=F,fit.file=NULL)

message("Fit the model for The Approx")
fit.approx <- fit.pw.parallel(init=init,datU = u.approx,coord=coord,reg=reg,reg.t=NULL,cutoff=cutoff,proppairs= 1,fixed=fixed,optim =T, hessian=F,sandwich=F,eps = 10^(-2), print.par.file=NULL,ncores=ncores,fit.load=F, fit.save=F,fit.file=NULL)

save.image(paste0("fit_simu_examples/fit_simu_",index,".RData"))












