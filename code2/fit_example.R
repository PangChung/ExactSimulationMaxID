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
n = 500
parR = c(1,0,10)
parGauss = c(1,10) 
x.coord <-  y.coord <- c(1:d)/(d+1) #grids 
coord = as.matrix(expand.grid(x.coord,y.coord))
Sigma = exp(-as.matrix(dist(coord))/1)*parR[3]

data.mh <- apply(mcmapply(mh,n=rep(1,n),SIMPLIFY = T,mc.cores = ncores,MoreArgs = list(parR=parR,Sigma=Sigma)),1,pG,parR=parR)
hist(data.mh[,100])
pairs <- t(combn(ncol(Sigma),2))

emp.pair(1000,pairs,0.9,-1/log(data.mh),parR=parR,cor.mat=Sigma/parR[3])



