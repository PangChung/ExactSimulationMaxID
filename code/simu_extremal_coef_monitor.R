### compute the high dimension extremal coefficients ###
### compute the empirical extremal coefficients ###
args <- commandArgs(trailingOnly=TRUE)
for(arg in args) {
  eval(parse(text = arg))
}
setwd("~/project/exact_simu_max_id_project/")
source("code_V3/Tools_MaxID.R");source("code_V3/Tools_Simu.R")

library(parallel)
library(evd)
library(mvtnorm)
library(ars)
#d = 7
D = d^2 # number of sites in the space
N = 100   # number of iterations when sample Y 
n = 10^4
LAMBDA.T = TRUE # whether the exponent measure is infinite or not on the whole space
parR = c(5,2) #alpha, beta
parGauss = c(2,-0.5,0,3)  #lambdas nu,
pars <- get.par(c(log(parR),parGauss),4)
x.coord <-  y.coord <- c(1:d)/(d+1) #grids 
cutoff <- 3/(7+1) # six order neighbors
coord = as.matrix(expand.grid(x.coord,y.coord))
coord = as.matrix(dist(coord))
reg=cbind(1,2*pnorm(coord[,1],0.5,0.25)-1)
reg.t = 0


file.list <- list.files(path="simu_extremal_coef/",pattern = paste0("simu_extremal_coef_",d,"_","\\d+.Rdata"),full.names = T)

while(length(file.list)<10){
  file.list <- list.files(path="simu_extremal_coef/",pattern = paste0("simu_extremal_coef_",d,"_","\\d+.Rdata"),full.names = T)
  Sys.sleep(60)
}

for(i in 1:length(file.list)){
  load(file.list[i],e<-new.env())
  if(is.list(e$u.ars)){
    ind <- sapply(e$u.ars,function(x){length(x)==d^2},simplify = T)
    e$u.ars <- matrix(unlist(e$u.ars[ind]),ncol=d^2,byrow = T)
  }
  if(is.list(e$u.mh)){
    ind <- sapply(e$u.mh,function(x){length(x)==d^2},simplify = T)
    e$u.mh <- matrix(unlist(e$u.mh[ind]),ncol=d^2,byrow = T)
  }
  if(is.list(e$u.approx)){
    ind <- sapply(e$u.approx,function(x){length(x)==d^2},simplify = T)
    e$u.approx <- matrix(unlist(e$u.approx[ind]),ncol=d^2,byrow = T)
  }
  save(list=ls(e),file=file.list[i],envir = e)
}

z <- qfrechet(c(0.05,0.25,0.5,0.75,0.95))
pairs2 <- t(combn(d^2,2))
pairs3 <- t(combn(d^2,3))

XDAT.list <- list(ars=NULL,mh=NULL,approx=NULL)
for(i in 1:length(file.list)){
  load(file.list[i],e<-new.env())
  XDAT.list$ars <- rbind(XDAT.list$ars,apply(e$u.ars,2,qfrechet,loc=0,scale=1,shape=1))
  XDAT.list$mh <- rbind(XDAT.list$mh,apply(e$u.mh,2,qfrechet,loc=0,scale=1,shape=1))
  XDAT.list$approx <- rbind(XDAT.list$approx,apply(e$u.approx,2,qfrechet,loc=0,scale=1,shape=1))
}

theta.dim2 <- theta.dim3 <- list()
for(i in 1:3){
  XDAT = XDAT.list[[i]]
  theta.dim2[[i]] <- mcmapply(emp.pair,k=1:nrow(pairs2),MoreArgs = list(pair=pairs2,XDAT=XDAT,pars=pars,reg=reg,reg.t=reg.t),SIMPLIFY = F,mc.cores = 40)
 # theta.dim3[[i]] <- mcmapply(emp.pair,k=1:nrow(pairs3),MoreArgs = list(pair=pairs3,XDAT=XDAT,pars=pars,reg=reg,reg.t=reg.t),SIMPLIFY = F,mc.cores = 40)
}

save.image(file=paste0("simu_extremal_coef/extremal_coef","_",d,".Rdata"))







