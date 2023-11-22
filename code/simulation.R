###control variables ### 
file.name="simulation.RData";original=T;N=10^5;index=0;year=1950
###load control variables ###
args <- commandArgs(trailingOnly=TRUE)
for(arg in args) {
  eval(parse(text = arg))
}
setwd("~/scratch/zhongp/exact_simu_max_id_project/code_V4")
source("Tools_MaxID.R");source("Tools_Simu.R")
library(parallel)
library(evd)
library(mgcv)
library(mvtnorm)
library(ars)
load("simulation.RData")
model.pred<-predict.gam(marginal.fit,data.frame(lon=data.info[stat2.avail,"LON"],lat=data.info[stat2.avail,"LAT"],t=rep(reg.t[year-1918+1],sum(stat2.avail)),alt=reg[stat2.avail,2]))
model.pred[,3]<-marginal.fit$fitted.values[1,3]
if(original){
  reg.t_model<-as.matrix(reg.t[year-1918+1])
  Z <- mcmapply(mh,n=1:N,SIMPLIFY = T,mc.cores = 64,MoreArgs = list(pars=par_model,reg=reg[stat2.avail,],reg.t=reg.t_model))
  Z <- mcmapply(function(i){pG(Z[,i],par_model$parR)},i=1:ncol(Z),SIMPLIFY = T,mc.cores = 64)
  Z  <- mcmapply(function(i){qgev(Z[,i],loc = model.pred[,1],scale = exp(model.pred[,2]),shape = model.pred[1,3])},i=1:ncol(Z),mc.cores=64)
  Z.summary = matrix(NA,nrow=ncol(Z),ncol=3)
  Z.summary[,1] <- mcmapply(function(id){min(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  Z.summary[,2] <-mcmapply(function(id){mean(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  Z.summary[,3]<-mcmapply(function(id){max(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  write.csv(Z.summary,file=paste0("/scratch/zhongp/project/output/simulation/simulation_results",year,"_",index,".RData"))
}else{
  par_model = get.par(fit.new$par,4)
  reg.t_model<-as.matrix(reg.t[year-1918+1])
  Z <- mcmapply(mh,n=1:N,SIMPLIFY = T,mc.cores = 64,MoreArgs = list(pars=par_model,reg=reg[stat2.avail,],reg.t=reg.t_model))
  Z <- mcmapply(function(i){pG(Z[,i],par_model$parR)},i=1:ncol(Z),SIMPLIFY = T,mc.cores = 64)
  Z  <- mcmapply(function(i){qgev(Z[,i],loc = model.pred[,1],scale = exp(model.pred[,2]),shape = model.pred[1,3])},i=1:ncol(Z),mc.cores=64)
  Z.summary = matrix(NA,nrow=ncol(Z),ncol=3)
  Z.summary[,1] <- mcmapply(function(id){min(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  Z.summary[,2] <-mcmapply(function(id){mean(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  Z.summary[,3]<-mcmapply(function(id){max(Z[,id])},id=1:ncol(Z),mc.cores=64,SIMPLIFY = T)
  write.csv(Z.summary,file=paste0("/scratch/zhongp/project/output/simulation/simulation_results",year,"_",index,".csv"),col.names = F,row.names = F)
}
