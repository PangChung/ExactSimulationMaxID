args <- commandArgs(trailingOnly=TRUE)
for(arg in args) {
  eval(parse(text = arg))
}
setwd("code2/")
source("Tools_MaxID.R");source("Tools_Simu.R")
library(parallel)
library(mvtnorm)
library(evd)
ncores=detectCores()
d = 10
D = d^2 # number of sites in the space
n = 100
parR = c(1,0,10)
x.coord <-  y.coord <- c(1:d)/(d+1) #grids 
coord = as.matrix(expand.grid(x.coord,y.coord))
Sigma = exp(-as.matrix(dist(coord))/0.5)*parR[3]

data.mh <- apply(mcmapply(mh,n=rep(1,n),SIMPLIFY = T,mc.cores = ncores,MoreArgs = list(parR=parR,Sigma=Sigma)),1,pG,parR=parR)
#hist(data.mh)
pairs <- t(combn(ncol(Sigma),2))
dist.vec = as.matrix(dist(coord))[pairs]

## extremal coefficient ##
## 2 is the index number of the pairs
## pairs contains all the pairs
## z=5 (default) is the level
## -1/log(data.mh) will transform the uniform data data.mh to the unit Frechet scale
## it produce a vector of length 3, (empirical estimates from the data, the standard deviation of the estimates, the theoretical true value)
emp.pair(99,pairs,z=0.5,-1/log(data.mh),parR=parR,cor.mat=Sigma/parR[3])
ext.est <- mcmapply(1:nrow(pairs),FUN=emp.pair,mc.cores=ncores,MoreArgs = list(z=1,pair=pairs,XDAT=-1/log(data.mh),parR=parR,cor.mat=Sigma/parR[3]),SIMPLIFY = TRUE)
plot(dist.vec,ext.est[1,],pch=20)
points(dist.vec,ext.est[3,],col=2,pch=20)

