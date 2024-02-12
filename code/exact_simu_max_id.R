setwd("code/")
source("Tools_MaxID.R");source("Tools_Simu.R")
library(parallel)
library(evd)
library(mvtnorm)
library(ars)
#### prepare the variables and data ####
d = 7
D = d^2 # number of sites in the space
N = 20  # number of iterations when sample Y 
mc.cores = 2
parR = c(5.5,2.4) #alpha, beta
parGauss = c(2.12,-0.31,0.23,3.2) #nu, lambdas
pars <- get.par(c(log(parR),parGauss),4)
x.coord <-  y.coord <- c(1:d)/(d+1) #grids 
cutoff <- 3/(7+1) # six order neighbors
coord = as.matrix(expand.grid(x.coord,y.coord))
coord = as.matrix(dist(coord))
reg=cbind(1,2*pnorm(coord[,1],0.5,0.25)-1)
reg.t = 0

t0 <- proc.time()
Y.approx <- mcmapply(rmaxidspat,n=rep(1,N),MoreArgs = list(coord=coord,parR=pars$parR,parGauss = pars$parGauss,reg=reg,reg.t=reg.t),SIMPLIFY = T,mc.cores = mc.cores)
print((t1 <- proc.time() - t0)/N*mc.cores)

t0 <- proc.time()
Y.mh <- mcmapply(mh,n=1:N,SIMPLIFY = T,mc.cores = mc.cores,MoreArgs = list(pars=pars,reg=reg,reg.t=reg.t))
print((t1 <- proc.time() - t0)/N*mc.cores)

t0 <- proc.time()
Y.ars <- mcmapply(mh,n=1:10^3,MoreArgs = list(ars=T,pars=pars,reg=reg,reg.t=reg.t),SIMPLIFY = T,mc.cores = mc.cores)
print((t1 <- proc.time() - t0)/N*mc.cores)

save.image()

load("simulation_2.RData",e<-new.env())
# ks.1<- sapply(1:49,function(n){ks.test(Y.ars[n,],Y.mh[n,])$p.value},simplify = T)
# ks.2<- sapply(1:49,function(n){ks.test(Y.ars[n,],Y.approx[n,])$p.value},simplify = T)
# ks.3<- sapply(1:49,function(n){ks.test(Y.approx[n,],Y.mh[n,])$p.value},simplify = T)
# ind <- sapply(e$Y.ars,length,simplify = T)==49
data.file <- list.files(path="simulation",full.names = T,pattern="simulation_")
ks.ars <- ks.mh <- ks.approx <- matrix(NA,nrow=49,ncol=length(data.file))
for(i in 1:length(data.file)){
load(data.file[i],e<-new.env())
#e$Y.ars <- matrix(unlist(e$Y.ars[ind]),nrow=49)
if(!is.list(e$Y.ars)){
e$u.ars <- apply(e$Y.ars,1,pG,parR=pars$parR)
e$u.mh <- apply(e$Y.mh,1,pG,parR=pars$parR)        
e$u.approx <- apply(e$Y.approx,1,pG,parR=pars$parR)
ks.ars[,i] <- apply(e$u.ars,2,function(x){ks.test(x,"punif")$p.value})
ks.mh[,i] <- apply(e$u.mh,2,function(x){ks.test(x,"punif")$p.value})
ks.approx[,i] <- apply(e$u.approx,2,function(x){ks.test(x,"punif")$p.value})
}
}

fun<-function(ks){
ind<-!is.na(ks)
return(sum(ks[ind] < 0.05)/sum(ind))
}

mean(apply(ks.ars,1,fun))
mean(apply(ks.mh,1,fun))
mean(apply(ks.approx,1,fun))


pdf(onefile = T,file="Simulation.pdf")
par(mfrow=c(1,1))
plot(1:(49*3),c(ks.ars,ks.mh,ks.approx),xaxt="n",ylim=c(0,1.5),col=rep(1:3,each=49),pch=rep(c(1:3),each=49),cex=0.5,
     xlab="site",ylab="ks-test p-value",main="KS test of the pseudo uniform samples against uniform distribution")
abline(h=0.05,lty=2,col="red")
abline(v=c(49.5,(49*2)+0.5),lty=2,col="red")
axis(1,at=c(1,49,49*2,49*3),labels = c("1","49 1","49 1","49"),tick = F)
legend("topright",col=1:3,pch=1:3,legend=c("ARS","MH","Original"),bty="n")

par(mfrow=c(3,1))
hist(c(u.ars),100,freq=F,xlab="",main="Pseudo uniform samples form Adaptive Rejection Sampling (ARS)")
abline(h=1,lty=2)
hist(c(u.mh),100,freq=F,xlab="",main = "Pseudo uniform samples form Metropolis-Hastings (MH)")
abline(h=1,lty=2)
hist(c(u.approx),100,freq=F,xlab="",main="Pseudo uniform samples form Original Sampling Method used in the paper")
abline(h=1,lty=2)
dev.off()

ks.ars.approx <- sapply(1:49,function(i){ks.test(u.ars[,i],u.approx[,i])$p.value})
ks.mh.approx <-  sapply(1:49,function(i){ks.test(u.mh[,i],u.approx[,i])$p.value})
ks.mh.ars <- sapply(1:49,function(i){ks.test(u.mh[,i],u.ars[,i])$p.value})
plot(c(ks.ars.approx,ks.mh.approx,ks.mh.ars),ylim=c(0,1),col=rep(1:3,each=49),pch=rep(c(1:3),each=49),cex=0.5)
abline(h=0.05)

par(mfrow=c(1,3))
image(x.coord,y.coord,matrix(Y.approx,nrow=d,ncol=d))
image(x.coord,y.coord,matrix(Y.mh,nrow=d,ncol=d))
image(x.coord,y.coord,matrix(Y.ars,nrow=d,ncol=d))

####Having Fun with ARS
# grad.new<-function(logends,x){
#   return(142/x-57)
# }
# #define the log concave funtion, which is the target log-form distribution.
# logends<-function(x){
#   return(142*log(x)-57*x)
# }
# nor2=my.ars(c(1,3),logends,0,Inf,1000)
# plot(density(nor2),xlab ="Values of Random Variables",ylab="Density",main = "Comparision of Sample and Its Distribution: Gamma Distribution")
# points(seq(min(nor2),max(nor2),0.001),dgamma(seq(min(nor2),max(nor2),0.001),141,57),pch=".",col="red")
# logends1<-function(x){
#   return(-1/2*x^2)
# }
# grad.new<-function(logends,x){
#   return(-x)
# }
# nor1=my.ars(c(-5,5),logends1,-Inf,Inf,1000)
# plot(density(nor1),xlab ="Values of Random Variables",ylab="Density",main = "Comparision of Sample and Its Distribution: Standard Guassian")
# points(seq(min(nor1),max(nor1),0.001),dnorm(seq(min(nor1),max(nor1),0.001)),pch=".",col="red")


### Fits on simulated data and plots ###
fit.files <- list.files(path="fit_simu_examples",pattern="fit_simu_*",full.names = T)
par.true = c(1,1,1,1,1)
par.ars <- par.mh <- par.approx <- matrix(NA,ncol=5,nrow=length(fit.files))
for(i in 1:length(fit.files)) {
load(fit.files[i],e<-new.env())
par.ars[i,] <- e$fit.ars$par[-5]
par.mh[i,] <- e$fit.mh$par[-5]
par.approx[i,] <- e$fit.approx$par[-5]
}
par.ars[,1:2] <- exp(par.ars[,1:2])
par.mh[,1:2] <- exp(par.mh[,1:2])
par.approx[,1:2] <- exp(par.approx[,1:2])
pdf(file="boxplot.pdf",onefile = T)
nn <- c("ARS","MH","Approx")
main = paste("Error w.r.t True value:",c("alpha","beta","lambda_0","lambda_1","nu"))
for(i in 1:5){
  boxplot(cbind(par.ars[,i],par.mh[,i],par.approx[,i])-par.true[i],names=nn,main=main[i])
  abline(h=0,col=2)
}
dev.off()

### plot the empirical extremal coef vs the theoretical one ###
#load("../data/extremal_coef_temp.Rdata")
load("../data/extremal_coef_temp.Rdata")
pdf("../pdf/scatter_plot_extremal.pdf",onefile = T,height=4,width = 4*3)
level <- c(0.05,0.25,0.5,0.75,0.95)
methd <- c("ARS","MH","APPROX")
## dim 2
for(i in 1:5){
  par(mfrow=c(1,3))
  for(j in 1:3){
    x <- sapply(theta.dim2[[j]],function(x){x[i,1]},simplify = T)
    y <- sapply(theta.dim2[[j]],function(x){x[i,3]},simplify = T)
    se <- sqrt(mean((x-y)^2))
    plot(x=x,y=y,cex=0.5,xlim=c(1,2),ylim=c(1,2),xlab="Empirical Extremal Coef",ylab="Theoretical Extremal Coef",main=paste0("Level: ",level[i]," Method: ",methd[j]," sd: ",round(se,4)))
    abline(a=0,b=1,col=2)
  }
}

## dim 3
for(i in 1:5){
  par(mfrow=c(1,3))
  for(j in 1:3){
    x <- sapply(theta.dim3[[j]],function(x){x[i,1]},simplify = T)
    y <- sapply(theta.dim3[[j]],function(x){x[i,3]},simplify = T)
    se <- sqrt(mean((x-y)^2))
    plot(x=x,y=y,cex=0.5,xlim=c(1,3),ylim=c(1,3),xlab="Empirical Extremal Coef",ylab="Theoretical Extremal Coef",main=paste0("Level: ",level[i]," Method: ",methd[j]," sd: ",round(se,4)))
    abline(a=0,b=1,col=2)
  }
}
dev.off()

########################################


