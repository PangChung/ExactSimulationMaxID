args <- commandArgs(trailingOnly=TRUE)
for(arg in args) {
  eval(parse(text = arg))
}

setwd("~/project/exact_simu_max_id_project/")
source("Tools_MaxID.R");source("Tools_Simu.R")
ncores = 40
library(parallel)
library(mvtnorm)
library(evd)

d = 7
D = d^2 # number of sites in the space
N = 100   # number of iterations when sample Y 
n = 10^3
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

Z.approx <- mcmapply(rmaxidspat,n=rep(1,n),MoreArgs = list(coord=coord,parR=pars$parR,parGauss = pars$parGauss,reg=reg,reg.t=reg.t),SIMPLIFY = F,mc.cores = ncores)

Z.mh <- mcmapply(mh,n=rep(1,n),SIMPLIFY = F,mc.cores = ncores,MoreArgs = list(reg=reg,reg.t=reg.t,pars=pars))

Z.ars <- mcmapply(mh,n=rep(1,n),MoreArgs = list(ars=T,reg=reg,reg.t=reg.t,pars=pars),SIMPLIFY = F,mc.cores = ncores)

save.image(file=paste0("simu_extremal_coef/simu_extremal_coef_",index,".Rdata"))

u.approx <- mcmapply(pG,x=Z.approx,MoreArgs = list(parR=pars$parR),SIMPLIFY = F,mc.cores = ncores)
u.mh <- mcmapply(pG,x=Z.mh,MoreArgs = list(parR=pars$parR),SIMPLIFY = F,mc.cores = ncores)
u.ars <- mcmapply(pG,x=Z.ars,MoreArgs = list(parR=pars$parR),SIMPLIFY = F,mc.cores = ncores)

save.image(file=paste0("simu_extremal_coef/simu_extremal_coef_",index,".Rdata"))


