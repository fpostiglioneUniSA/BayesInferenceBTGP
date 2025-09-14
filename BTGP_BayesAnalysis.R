library("R2OpenBUGS")
library(coda)
if(Sys.info()["sysname"] == "Windows") require(devEMF)

# For reproducibility purposes
set.seed(2023)

# data of liners
dataLiners <- list(
tempi1 = c(0, 11300, 14680, 31270), 
tempi2 = c(0, 11300, 21970), 
tempi3 = c(0, 12300, 16300), 
tempi4 = c(0, 14810, 18700, 28000), 
tempi5 = c(0, 10000, 30450, 37310), 
tempi6 = c(0, 6860, 17200, 24710), 
tempi7 = c(0, 2040, 12580, 16620), 
tempi8 = c(0, 7540, 8840, 9770, 16300), 
lin1 = c(0, 0.9, 1.3, 2.85),
lin2 = c(0, 1.5, 2),
lin3 = c(0, 1, 1.35),
lin4 = c(0, 1.9, 2.25, 2.75),
lin5 = c(0, 1.2, 2.75, 3.05),
lin6 = c(0, 0.5, 1.45, 2.15),
lin7 = c(0, 0.4, 2, 2.35),
lin8 = c(0, 0.5, 1.1, 1.15, 2.1)
)

DD = 4.0   # Maximum wear

numlin=8   # liners number

# Modify the name of the .RData containing the MCMC sample to analyze
#load("xxx.RData")



# --------------------------------
# posterior PDFs of parameters 
# --------------------------------
if (codpkg==FALSE) {
  if (nochains==3) {
    mc.sel <- rbind(mc.out1,mc.out2,mc.out3)
  } else if (nochains==2) {
    mc.sel <- rbind(mc.out1,mc.out2)
  } else mc.sel <- mc.out1
} else {
  if (nochains==3) {
    mc.sel <- rbind(mc.out1[[1]],mc.out2[[1]],mc.out3[[1]])
  } else if (nochains==2) {
    mc.sel <- rbind(mc.out1[[1]],mc.out2[[1]])
  } else mc.sel <- mc.out1[[1]]
}
simnum <- nrow(mc.sel) # consider the proper data size after combining MCMC chains 

cat("\n Data size after combining MCMC chains:",simnum,"\n")
print(summary(mc.sel))

x11()
#dev.new()
h<-hist(mc.sel[,"b"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h$breaks[-length(h$breaks)]+h$breaks[-1]),h$density, main="Pdf of b", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("bpostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h2<-hist(mc.sel[,"wlim"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h2$breaks[-length(h2$breaks)]+h2$breaks[-1]),h2$density, main="Pdf of wlim", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("wlimpostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h3<-hist(mc.sel[,"a"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h3$breaks[-length(h3$breaks)]+h3$breaks[-1]),h3$density, main="Pdf of a", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("apostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h4<-hist(mc.sel[,"beta"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h4$breaks[-length(h4$breaks)]+h4$breaks[-1]),h4$density, main="Pdf of beta", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("betapostpdf_",suffnome,".pdf"))

# # Export in csv (to be modified to be converted into .xls)
# write.csv(cbind(1/2*(h$breaks[-length(h$breaks)]+h$breaks[-1]),h$density), file = "bpost.csv",row.names = F)
# write.csv(cbind(1/2*(h2$breaks[-length(h2$breaks)]+h2$breaks[-1]),h2$density), file = "wlim_Upost.csv",row.names = F)
# write.csv(cbind(1/2*(h3$breaks[-length(h3$breaks)]+h3$breaks[-1]),h3$density), file = "apost.csv",row.names = F)
# write.csv(cbind(1/2*(h4$breaks[-length(h4$breaks)]+h4$breaks[-1]),h4$density), file = "betapost.csv",row.names = F)

# Traceplots (and density plots)
# install.packages("basicMCMCplots")
# samplelist<- list(codaobject1,codaobject2,codaobject3)
# samplelist<- list(codaobject1,codaobject2) # 2 chains
samplelist<- codaobject1 # 1 chain 

library(basicMCMCplots)
samplesPlot(samplelist)
dev.print(device=pdf, "Trace_densityPlots.pdf")
dev.print(device=png, "Trace_densityPlots.png", width=1000,height=1000)
dev.off()
chainsPlot(samplelist)
dev.new()
chainsPlot(samplelist,densityplot = F,height=1000,width=1000)


#

cat("Mean of b: ", mean(mc.sel[,'b']), "\n")
cat("Mean of wlim = U", mean(mc.sel[,'wlim']),"\n")
cat("Mean of a", mean(mc.sel[,'a']),"\n")
cat("Mean of beta", mean(mc.sel[,'beta']),"\n")




# HPD intervals

print(" ")
print("0.90 HPD interval for b")
print(" ")
Hpdb <- HPDinterval(as.mcmc(mc.sel[,'b']),  prob = 0.9)
print(Hpdb[1,])
print(" ")

print(" ")
print("0.90 HPD interval for wlim")
print(" ")
Hpdwlim <- HPDinterval(as.mcmc(mc.sel[,'wlim']),  prob = 0.9)
print(Hpdwlim[1,])
print(" ")

print(" ")
print("0.90 HPD interval for a")
print(" ")
Hpda <- HPDinterval(as.mcmc(mc.sel[,'a']),  prob = 0.9)
print(Hpda[1,])
print(" ")

print(" ")
print("0.90 HPD interval for beta")
print(" ")
Hpdbet <- HPDinterval(as.mcmc(mc.sel[,'beta']),  prob = 0.9)
print(Hpdbet[1,])
print(" ")




# ------------------------------------
# W(t): mean and variance
# ------------------------------------

cat("\n #-------------------------\n")
cat(" # W(t): mean and variance \n")
cat(" #-------------------------\n")



Ew <- function(t, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  integrand <- function(w){
    ifelse(w<wli,w*(bet*wli)/(wli-w)^2*((bet*w)/(wli-w))^((t/aa)^bb-1)/gamma((t/aa)^bb)*exp(-(bet*w)/(wli-w)),0)
    }
  resul <- integrate(integrand, lower=0, upper = wli, stop.on.error=F)$value
  return(resul)
}

Ew2 <- function(t, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  integrand <- function(w){
    ifelse(w<wli,1-pgamma(bet*w/(wli-w), shape = (t/aa)^bb),0)
  }
  resul <- integrate(integrand, lower=0, upper = Inf, stop.on.error=F)$value
  return(resul)
}

EwSquared <- function(t, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  integrand <- function(w){
    ifelse(w<wli,w*(1-pgamma(bet*w/(wli-w), shape = (t/aa)^bb)),0)
  }
  resul <- 2*integrate(integrand, lower=0, upper = Inf, stop.on.error=F)$value
  return(resul)
}


# initialize matrices for mean W(t)
et <- seq(10,80000,1000)
matr_ew <- matrix(,nrow = nrow(mc.sel), ncol = length(et))
ew <- rep(0,length(et))
ewhpd <- matrix(,nrow=length(et), ncol=2)
matr_ew2 <- matrix(,nrow = nrow(mc.sel), ncol = length(et))
ew2 <- rep(0,length(et))
ewhpd2 <- matrix(,nrow=length(et), ncol=2)


# initialize matrices for var W(t)
vt <- seq(10,80000,1000)
matr_var <- matrix(,nrow = nrow(mc.sel), ncol = length(vt))
vw <- rep(0,length(vt))
vwhpd <- matrix(,nrow=length(vt), ncol=2)


for(jj in 1:length(et)) {
  for(ii in 1:nrow(mc.sel)){matr_ew2[ii,jj]<-Ew2(et[jj],wli=mc.sel[ii,"wlim"], bet=mc.sel[ii,"beta"], aa=mc.sel[ii,"a"], bb=mc.sel[ii,"b"])}
  ewhpd2[jj,] <- HPDinterval(as.mcmc(matr_ew2[,jj]), prob=0.9)[1,]
}
ew2<-colMeans(matr_ew2); print(ew2)


for(jj in 1:length(vt)) {
  for(ii in 1:nrow(mc.sel)){matr_var[ii,jj]<-EwSquared(vt[jj],wli=mc.sel[ii,"wlim"], bet=mc.sel[ii,"beta"], aa=mc.sel[ii,"a"], bb=mc.sel[ii,"b"]) - matr_ew2[ii,jj]^2}
  vwhpd[jj,] <- HPDinterval(as.mcmc(matr_var[,jj]), prob=0.9)[1,]
}
varw<-colMeans(matr_var); print(varw)


# Empirical mean & variance of wear
EmpMean <- function(tau) {
  indexSum <- rep(0,numlin)
  for(ii in 1:numlin){indexSum[ii]=max(dataLiners[[ii]])>=tau}
  ii=1
  result <- 0
  while(ii<=numlin) {
    if(indexSum[ii]!=0){
      i2 <- min(which(dataLiners[[ii]]>=tau))
      i1 <- i2 - 1
      xx <- c(dataLiners[[ii]][i1],dataLiners[[ii]][i2])
      yy <- c(dataLiners[[ii+numlin]][i1],dataLiners[[ii+numlin]][i2])
      result <- result + approx(xx,yy,xout=tau)$y
    }
    ii=ii+1
  }
  return(result/sum(indexSum))
}

EmpVar <- function(tau) {
  indexSum <- rep(0,numlin)
  for(ii in 1:numlin){indexSum[ii]=max(dataLiners[[ii]])>=tau}
  ii=1
  result <- 0
  while(ii<=numlin) {
    if(indexSum[ii]!=0){
      i2 <- min(which(dataLiners[[ii]]>=tau))
      i1 <- i2 - 1
      xx <- c(dataLiners[[ii]][i1],dataLiners[[ii]][i2])
      yy <- c(dataLiners[[ii+numlin]][i1],dataLiners[[ii+numlin]][i2])
      result <- result + (approx(xx,yy,xout=tau)$y - EmpMean(tau))^2
    }
    ii=ii+1
  }
  return(result/(sum(indexSum)-1))
}

taus<-seq(2500,35000,by=2500)
emean = rep(0,length(taus))
evar = rep(0,length(taus))
for(ii in 1:length(taus)){
  emean[ii] <- EmpMean(taus[ii]) 
  evar[ii] <- EmpVar(taus[ii]) 
}

taus <- c(0,taus)
emean <- c(0,emean)
evar <- c(0,evar)
# Export in csv (to be modified to be converted into .xls)
write.csv(cbind(taus,emean), file = "empMean_w.csv", row.names=FALSE)
write.csv(cbind(taus,evar), file = "empVar_w.csv", row.names=FALSE)

x11() 
#dev.new() # Mean wear and 0.90 HPD
plot(et,ew2, main="Mean wear [mm]", xlab="Operating time t [h]", ylab="Mean wear [mm]", xlim=c(0,80000), ylim=c(0,5.1), type="l", col="red",lwd=2)
matplot(et,ewhpd2[,2], pch=2, lty=2, lwd=1, add=T, col="red", type="l")
matplot(et,ewhpd2[,1], pch=2, lty=2, lwd=1, add=T, col="red", type="l")
matplot(taus,emean, pch=16, add=T, col="black", type="p")
legend('bottomright', inset = 0.05, c('Bayesian estimate',"0.9 HPD limits",NA,"empirical estimate"), lty=c(1,2,NA,NA), lwd=c(2,1,1,NA), col=c('red','red',NA, 'black'), pch=c(NA,NA,NA,16), cex=1)
dev.print(device=pdf, "MeanWear.pdf")

x11()
#dev.new() # Mean wear and 0.90 HPD
plot(vt,varw, main="Variance of wear [mm^2]", xlab="Operating time t [h]", ylab="Variance of wear [mm^2]", xlim=c(0,80000), ylim=c(0,0.16), type="l", col="red",lwd=2)
matplot(vt,vwhpd[,2], pch=2, lty=2, lwd=1, add=T, col="red", type="l")
matplot(vt,vwhpd[,1], pch=2, lty=2, lwd=1, add=T, col="red", type="l")
matplot(taus,evar, pch=16, add=T, col="black", type="p")
legend('topright', inset = 0.05, c('Bayesian estimate',"0.9 HPD limits",NA,"empirical estimate"), lty=c(1,2,NA,NA), lwd=c(2,1,1,NA), col=c('red','red',NA, 'black'), pch=c(NA,NA,NA,16), cex=1)
dev.print(device=pdf, "VarianceWear.pdf")



#----------------------------------------------------
# Reliability function R(x) and cond. rel. R_t(x|w_t)
#----------------------------------------------------

cat(sprintf("\n #-------------------------------------\n"))
cat(sprintf(" # Reliability and cond. rel. functions\n"))
cat(sprintf(" #-------------------------------------\n"))


Rxw <- function(x, t, wt, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  resul <- pgamma(bet*(wm/(wli-wm)-wt/(wli-wt)), shape = ((t+x)/aa)^bb - (t/aa)^bb)
  return(resul)
}

RelFun <- function(t, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  resul <- pgamma(bet*wm/(wli-wm), shape = (t/aa)^bb)
  return(resul)
} # from Rxw definition: resul <- Rxw(x,t=0,wt=0)

RelFunGEL <- function(t, pGEL=1, wm=DD, wli=mc.sel[,"wlim"],  bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  resul <- mean(RelFun(t,wm,wli,bet,aa,bb)^(-pGEL))^(-1/pGEL)
  return(resul)
}

xmax <- 200000
xstep <- xmax/50
xxR <- seq(0,xmax,xstep)
Rfun <- rep(0,length(xxR))
RfunGEL <- rep(0,length(xxR))
Rfunhpd <- matrix(,nrow=length(xxR), ncol=2)
Rfunci <- matrix(,nrow=length(xxR), ncol=2)
Rfunlow <- rep(0,length(xxR))

GELratio <- function(d,delta,rr){return((delta^d - d*log(delta) - 1)/(delta^(-d) + d*log(delta) - 1) - rr)}

# pGELrel <- 2.87 # GEL parameter r = 1.2
# pGELrel <- 1.50 # GEL parameter r = 1.1
pGELrel <- 0.784 # = round(resRel$root,digits=2) # GEL parameter r = 1.1 e delta = 1.2

Sys.time()->tstart
for(ii in 1:length(xxR)) {
  Rfun[ii] <- mean(RelFun(xxR[ii], wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]))
  RfunGEL[ii] <- RelFunGEL(xxR[ii],pGEL=pGELrel)
}


x11()
#dev.new() 
#plot(xxR,RfunGEL, main="Reliability", xlab="time", ylab="Reliability", xlim=c(0,xmax), ylim=c(0,1.05), type="l", col="black")
plot(xxR,Rfun, main="Reliability", xlab="time", ylab="Reliability", xlim=c(0,xmax), ylim=c(0,1.05), type="l", col="black")
matplot(xxR,RfunGEL, pch=2, lty=3, add=T, col="blue", type="l")
legend('bottomleft', c('R(t) SEL est.','R(t) GEL est.'), lty=c(1,3), col=c('black','blue'), cex=.9)
dev.print(device=pdf, "RelFun.pdf")

#-------------------------------------
#  Conditional residual reliability
#-------------------------------------

cat(sprintf("\n #-------------------------------------------\n"))
cat(sprintf(" # Conditional residual reliability R_t(x|wt) \n"))
cat(sprintf(" #-------------------------------------------\n"))

Rxw <- function(x, t, wt, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  resul <- pgamma(bet*(wm/(wli-wm)-wt/(wli-wt)), shape = ((t+x)/aa)^bb - (t/aa)^bb)
  return(resul)
} # like above

RxwGEL <- function(x,t, wt, pGEL=1, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  resul <- mean(Rxw(x,t,wt,wm,wli,bet,aa,bb)^(-pGEL))^(-1/pGEL)
  return(resul)
}

xmax <- 200000; xstep <- xmax/50
xxR <- seq(0,xmax,xstep)
Rx1 <- rep(0,length(xxR))
Rx2 <- rep(0,length(xxR))
Rx3 <- rep(0,length(xxR))
Rx4 <- rep(0,length(xxR))
Rx5 <- rep(0,length(xxR))
Rx6 <- rep(0,length(xxR))
Rx7 <- rep(0,length(xxR))
Rx8 <- rep(0,length(xxR))
Rx1GEL <- rep(0,length(xxR))
Rx2GEL <- rep(0,length(xxR))
Rx3GEL <- rep(0,length(xxR))
Rx4GEL <- rep(0,length(xxR))
Rx5GEL <- rep(0,length(xxR))
Rx6GEL <- rep(0,length(xxR))
Rx7GEL <- rep(0,length(xxR))
Rx8GEL <- rep(0,length(xxR))
Rx1hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx2hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx3hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx4hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx5hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx6hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx7hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx8hpd <- matrix(,nrow=length(xxR), ncol=2)
Rx1ci <- matrix(,nrow=length(xxR), ncol=2)
Rx2ci <- matrix(,nrow=length(xxR), ncol=2)
Rx3ci <- matrix(,nrow=length(xxR), ncol=2)
Rx4ci <- matrix(,nrow=length(xxR), ncol=2)
Rx5ci <- matrix(,nrow=length(xxR), ncol=2)
Rx6ci <- matrix(,nrow=length(xxR), ncol=2)
Rx7ci <- matrix(,nrow=length(xxR), ncol=2)
Rx8ci <- matrix(,nrow=length(xxR), ncol=2)


for(ii in 1:length(xxR)) {
  # Rx1[ii] <- mean(Rxw(xxR[ii],dataLiners[[1]][4],dataLiners[[9]][4]))
  # Rx2[ii] <- mean(Rxw(xxR[ii],dataLiners[[2]][3],dataLiners[[10]][3]))
  Rx3[ii] <- mean(Rxw(xxR[ii],dataLiners[[3]][3],dataLiners[[11]][3]))
  # Rx4[ii] <- mean(Rxw(xxR[ii],dataLiners[[4]][4],dataLiners[[12]][4]))
  Rx5[ii] <- mean(Rxw(xxR[ii],dataLiners[[5]][4],dataLiners[[13]][4]))
  # Rx6[ii] <- mean(Rxw(xxR[ii],dataLiners[[6]][4],dataLiners[[14]][4]))
  # Rx7[ii] <- mean(Rxw(xxR[ii],dataLiners[[7]][4],dataLiners[[15]][4]))
  # Rx8[ii] <- mean(Rxw(xxR[ii],dataLiners[[8]][5],dataLiners[[16]][5]))
  # Rx1GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[1]][4],dataLiners[[9]][4],pGEL=pGELrel)
  # Rx2GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[2]][3],dataLiners[[10]][3],pGEL=pGELrel)
  Rx3GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[3]][3],dataLiners[[11]][3],pGEL=pGELrel)
  # Rx4GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[4]][4],dataLiners[[12]][4],pGEL=pGELrel)
  Rx5GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[5]][4],dataLiners[[13]][4],pGEL=pGELrel)
  # Rx6GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[6]][4],dataLiners[[14]][4],pGEL=pGELrel)
  # Rx7GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[7]][4],dataLiners[[15]][4],pGEL=pGELrel)
  # Rx8GEL[ii] <- RxwGEL(xxR[ii],dataLiners[[8]][5],dataLiners[[16]][5],pGEL=pGELrel)
}



# Export in csv (to be modified to be converted into .xls)
write.csv(cbind(xxR,Rx3,Rx3GEL), file = "R3.csv", row.names=FALSE)
write.csv(cbind(xxR,Rx5,Rx5GEL), file = "R5.csv", row.names=FALSE)

x11()
plot(xxR,Rx3, main="Conditional residual reliability R_t(x|wt)", xlab="RUL x", ylab="R_{t}(x|w_t)", xlim=c(0,xmax), ylim=c(0,1.05), type="l", col=1)
matplot(xxR,Rx5, pch=2, lty=1, add=T, col=5, type="l")
matplot(xxR,Rx3GEL, pch=2, lty=3, add=T, col=1, type="l")
matplot(xxR,Rx5GEL, pch=2, lty=3, add=T, col=5, type="l")
legend('topright', c('RUL SEL 3','RUL SEL 5','RUL GEL 3','RUL GEL 5'), lty=c(1, 1, 3, 3), col=c(1,5,1,5), cex=.9)
dev.print(device=pdf, "CondResRel.pdf")



# ---------------------------------------
# Degradation increment pdf DeltaW
# ---------------------------------------

cat(sprintf("\n #----------------------------------\n"))
cat(sprintf(" # Degradation increment DeltaW pdf \n"))
cat(sprintf(" #----------------------------------\n"))

ww <- seq(0,2,0.025)
w1<-w2<-w3<-w4<-w5<-w6<-w7<-w8 <- rep(0,length(ww))


# Sampling from posterior distrib. of the wear incrememnt dW
dW <- function(t, wt, dT = 1E4, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
  zwt <- rgamma(simnum,shape = ((t + dT)/aa)^bb - (t/aa)^bb, scale = 1)
  resul <- wli/(bet/(zwt+bet*wt/(wli-wt)) + 1) - wt
  return(resul)
}

Sys.time()->tstart
dWnew <- dW(0,0)
dW1 <- dW(dataLiners[[1]][4],dataLiners[[9]][4])
dW2 <- dW(dataLiners[[2]][3],dataLiners[[10]][3])
dW3 <- dW(dataLiners[[3]][3],dataLiners[[11]][3])
dW4 <- dW(dataLiners[[4]][4],dataLiners[[12]][4])
dW5 <- dW(dataLiners[[5]][4],dataLiners[[13]][4])
dW6 <- dW(dataLiners[[6]][4],dataLiners[[14]][4])
dW7 <- dW(dataLiners[[7]][4],dataLiners[[15]][4])
dW8 <- dW(dataLiners[[8]][5],dataLiners[[16]][5])


# dW pdfs
hWnew <- hist(dWnew, breaks = c(seq(0,2.5,.05),max(dWnew)+.1), plot=F)
hW1 <- hist(dW1, breaks = c(seq(0,1.25,.025),max(dW1)+.1), plot=F)
hW2 <- hist(dW2, breaks = c(seq(0,25,.05),max(dW2)+.1), plot=F)
hW3 <- hist(dW3, breaks = c(seq(0,25,.05),max(dW3)+.1), plot=F)
hW4 <- hist(dW4, breaks = c(seq(0,1.25,.025),max(dW4)+.1), plot=F)
hW5 <- hist(dW5, breaks = c(seq(0,1.25,.025),max(dW5)+.1), plot=F)
hW6 <- hist(dW6, breaks = c(seq(0,25,.05),max(dW6)+.1), plot=F)
hW7 <- hist(dW7, breaks = c(seq(0,25,.05),max(dW7)+.1), plot=F)
hW8 <- hist(dW8, breaks = c(seq(0,25,.05),max(dW8)+.1), plot=F)

Edwnew <- mean(dWnew)
Edw1 <- mean(dW1)
Edw2 <- mean(dW2)
Edw3 <- mean(dW3)
Edw4 <- mean(dW4)
Edw5 <- mean(dW5)
Edw6 <- mean(dW6)
Edw7 <- mean(dW7)
Edw8 <- mean(dW8)

cat("Mean dW after 10000h for new unit:",Edwnew,"\n")
cat("Mean dW after 10000h for liner 1:",Edw1,"\n")
cat("Mean dW after 10000h for liner 2:",Edw2,"\n")
cat("Mean dW after 10000h for liner 3:",Edw3,"\n")
cat("Mean dW after 10000h for liner 4:",Edw4,"\n")
cat("Mean dW after 10000h for liner 5:",Edw5,"\n")
cat("Mean dW after 10000h for liner 6:",Edw6,"\n")
cat("Mean dW after 10000h for liner 7:",Edw7,"\n")
cat("Mean dW after 10000h for liner 8:",Edw8,"\n\n")

# 0.9 HPD intervals for dWs
HpddWnew <- HPDinterval(as.mcmc(dWnew),  prob = 0.9)
HpddW1 <- HPDinterval(as.mcmc(dW1),  prob = 0.9)
HpddW2 <- HPDinterval(as.mcmc(dW2),  prob = 0.9)
HpddW3 <- HPDinterval(as.mcmc(dW3),  prob = 0.9)
HpddW4 <- HPDinterval(as.mcmc(dW4),  prob = 0.9)
HpddW5 <- HPDinterval(as.mcmc(dW5),  prob = 0.9)
HpddW6 <- HPDinterval(as.mcmc(dW6),  prob = 0.9)
HpddW7 <- HPDinterval(as.mcmc(dW7),  prob = 0.9)
HpddW8 <- HPDinterval(as.mcmc(dW8),  prob = 0.9)
print(" ")
print("0.9 HPD interval after 10000h for a new unit")
print(HpddWnew[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 1")
print(HpddW1[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 2")
print(HpddW2[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 3")
print(HpddW3[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 4")
print(HpddW4[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 5")
print(HpddW5[1,])
print("0.9 HPD interval after 10000h for liner 6")
print(HpddW6[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 7")
print(HpddW7[1,])
print(" ")
print("0.9 HPD interval after 10000h for liner 8")
print(HpddW8[1,])
print(" ")

GELratio <- function(d,delta,rr){return((delta^d - d*log(delta) - 1)/(delta^(-d) + d*log(delta) - 1) - rr)}
res <- uniroot(GELratio, interval = c(-10, 10), delta=1.1,rr=0.8)

GELest <- function(x,pGEL=1) {
  resul <- mean(x^(-pGEL))^(-1/pGEL)
  return(resul)
}

# pGELlife <- 2.16 # GEL parameter for mean lifetime and RUL
# pGELlife <- 2.87 # GEL parameter for mean lifetime and RUL
pGELdW <- round(res$root,digits=2) # GEL parameter for wear increment, dW = -3.51 for r=0.8 & delta=1.1

# GEL dW
dWGELnew <- GELest(dWnew,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW of a new unit by generating a sample: %f h\n\n",dWGELnew))
dWGEL1 <- GELest(dW1,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW1 by generating a sample: %f h\n\n",dWGEL1))
dWGEL2 <- GELest(dW2,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW2 by generating a sample: %f h\n\n",dWGEL2))
dWGEL3 <- GELest(dW3,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW3 by generating a sample: %f h\n\n",dWGEL3))
dWGEL4 <- GELest(dW4,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW4 by generating a sample: %f h\n\n",dWGEL4))
dWGEL5 <- GELest(dW5,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW5 by generating a sample: %f h\n\n",dWGEL5))
dWGEL6 <- GELest(dW6,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW6 by generating a sample: %f h\n\n",dWGEL6))
dWGEL7 <- GELest(dW7,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW7 by generating a sample: %f h\n\n",dWGEL7))
dWGEL8 <- GELest(dW8,pGEL=pGELdW)
cat(sprintf("\nGEL est. of dW8 by generating a sample: %f h\n\n",dWGEL8))

write.csv(cbind(hWnew$mids,hWnew$density), file = "deltaW_new.csv", row.names=FALSE)
write.csv(cbind(hW3$mids,hW3$density), file = "deltaW_lin3.csv", row.names=FALSE)
write.csv(cbind(hW5$mids,hW5$density), file = "deltaW_lin5.csv", row.names=FALSE)
# write.csv(cbind(hW8$mids,hW8$density), file = "hW8.csv", row.names=FALSE)

x11()
#dev.new() 
plot(hWnew$mids,hWnew$density, pch=2, lty=2, main="Degradation increment pdf", xlab="w", ylab="f_DW(w)", xlim=c(0,2.5), ylim=c(0,6.0), type="l", col=1)
# matplot(hW1$mids,hW1$density, pch=3, lty=1, add=T, type="l", col=1)
# matplot(hW2$mids,hW2$density, pch=3, lty=1, add=T, type="l", col=2)
matplot(hW3$mids,hW3$density, pch=3, lty=1, add=T, type="l", col=3)
# matplot(hW4$mids,hW4$density, pch=3, lty=1, add=T, type="l", col=4)
matplot(c(hW5$mids,2.45),c(hW5$density,0), pch=3, lty=1, add=T, type="l", col=5)
# matplot(hW6$mids,hW6$density, pch=3, lty=1, add=T, type="l", col=6)
# matplot(hW7$mids,hW7$density, pch=3, lty=1, add=T, type="l", col=7)
# matplot(hW8$mids,hW8$density, pch=4, lty=1, add=T, type="l", col=8)
# legend('topright', c('new unit','item 1','item 2','item 3','item 4','item 5','item 6','item 7','item 8'), lty=c(2,rep(1,8)), col=c(1,1:8), cex=.9)
legend('topright', c('new unit','item 3','item 5'), lty=c(2,rep(1,2)), col=c(1,3,5), cex=.9)
dev.print(device=pdf, "DegrIncr.pdf")



# ------------------------------------
# Lifetime and RUL prediction
# ------------------------------------

cat(sprintf("\n #-----------------------------------------\n"))
cat(" # Lifetime X and RUL Xt|wt prediction pdf \n")
cat(sprintf(" #-----------------------------------------\n"))

FX <- function(x, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
 resul <- pgamma(bet*wm/(wli-wm), shape = (x/aa)^bb)
 return(resul)
} # CCDF of X <-> R(x) with g(w)=beta*w/(wlim-w)


FXtw <- function(x, t, wt, wm=DD, wli=mc.sel[,"wlim"], bet=mc.sel[,"beta"], aa=mc.sel[,"a"], bb=mc.sel[,"b"]) {
 resul <- pgamma(bet*(wm/(wli-wm)-wt/(wli-wt)), shape = ((t+x)/aa)^bb - (t/aa)^bb)
  return(resul)
} # CCDF of X_t|wt <-> Rt(x|wt) with g(w)=beta*w/(wlim-w)

GELest <- function(x,pGEL=1) {
  resul <- mean(x^(-pGEL))^(-1/pGEL)
  return(resul)
}


# pGELlife <- 2.87 # old GEL parameter for mean lifetime and RUL
pGELlife <- 2.16 # = round(resLife$root,digits=2) # GEL parameter for mean lifetime and RUL (\delta=1.2 and \rho=1.3)

Xm <- numeric(simnum)
Xmtw <- matrix(0,nrow=simnum,ncol=numlin) # numlin = 8 items
# Compute the mean of lifetime X for each posterior vector element
meanX <- function(wm, wli, bet, aa, bb) {
  integrand <- function(x){
    RelFun(x,wm,wli,bet,aa,bb)
  }
  resul <- integrate(integrand, lower=0, upper = 1e7, stop.on.error=F)$value
  return(resul)
}

for(ii in 1:simnum) Xm[ii]<-meanX(wm=DD,wli=mc.sel[ii,"wlim"], bet=mc.sel[ii,"beta"], aa=mc.sel[ii,"a"], bb=mc.sel[ii,"b"])


# CAMPIONAMENTO DA FX
FX.inv <- function(y, wm, wli, bet, aa, bb){uniroot(function(x){FX(x, wm, wli, bet, aa, bb) - y},interval=c(0,3000000))$root}
FX.inv <- Vectorize(FX.inv)
Sys.time() -> tstart
Un <- runif(simnum,0,1)   # random sample from U[0,1]
XX <- FX.inv(Un, DD, mc.sel[1:length(Un),"wlim"], mc.sel[1:length(Un),"beta"], mc.sel[1:length(Un),"a"], mc.sel[1:length(Un),"b"])
tTCDF <- (Sys.time()-tstart)
cat("\n Tempo generazione campione lifetime XX tramite inv CDF:",tTCDF, "s.\n")

h.XX <- hist(XX, breaks = c(seq(0,200000,2000),max(XX)+1), plot=F)

# 3. SEL X
XXSEL <- mean(XX) # stima SEL vita X == media X
cat(sprintf("\nSEL est. of life X by generating a sample: %f h\n\n",XXSEL))
HpdX <- HPDinterval(as.mcmc(XX),  prob = 0.9)
cat("0.9 HPD interval of lifetime X\n")
print(HpdX[1,])

# 3. GEL X
XXGEL <- GELest(XX,pGEL=pGELlife)
cat(sprintf("\nGEL est. of life X by generating a sample: %f h\n\n",XXGEL))


# CAMPIONAMENTO DA FXtw
XXtw <- matrix(0,nrow=simnum,ncol=numlin) # numlin = 8 items
FXtw.inv <- function(y, t, wt, wm, wli, bet, aa, bb){uniroot(function(x){FXtw(x, t, wt, wm, wli, bet, aa, bb) - y},interval=c(0,30000000))$root}
FXtw.inv <- Vectorize(FXtw.inv)
Sys.time() -> tstart
for(ii in 1:numlin){
  Untw <- runif(simnum,0,1)   # random sample from U[0,1]
  XXtw[,ii] <- FXtw.inv(Untw, dataLiners[[ii]][length(dataLiners[[ii]])],dataLiners[[ii+numlin]][length(dataLiners[[ii]])], DD, mc.sel[1:length(Untw),"wlim"], mc.sel[1:length(Untw),"beta"], mc.sel[1:length(Untw),"a"], mc.sel[1:length(Untw),"b"])
}

XXtwGEL <- integer(numlin)
XXtwSEL <- integer(numlin)
HpdXXtw <- matrix(0,nrow=numlin,ncol=2)
h.XXtw <- list()

for(ii in 1:numlin){
  XXtwSEL[ii] <- mean(XXtw[,ii])  # stima sulla v.a. predetta di Xtw 
  XXtwGEL[ii] <- GELest(XXtw[,ii],pGEL=pGELlife) # stima GEL su Xtw 
  # ETtw[ii] <- mean(Tmtw[,ii])  # stima sulla v.a. predetta vita residua
  # TTGELtw[ii] <- GELest(TTtw[,ii],pGEL=pGELlife)
  cat(sprintf("SEL est. of RUL X|wt per item #%d: %f\n",ii,XXtwSEL[ii]))
  cat(sprintf("GEL est. (p=%f) of RUL X|wt per item #%d: %f\n",pGELlife,ii,XXtwGEL[ii]))
  HpdXXtw[ii,] <- HPDinterval(as.mcmc(XXtw[,ii]),  prob = 0.9)[1,]
  cat(sprintf("0.9 HPD interval of RUL X|wt per item #%d:\n",ii))
  print(HpdXXtw[ii,])
  cat("\n")
  h.XXtw <- append(h.XXtw, list(hist(XXtw[,ii], breaks = c(seq(0,200000,2000),max(XXtw[,ii])+1), plot=F)))
}


# Export in csv (to be modified to be converted into .xls)
write.csv(cbind(1/2*(h.XX$breaks[-length(h.XX$breaks)]+h.XX$breaks[-1]),h.XX$density), file = "Xpdf.csv", row.names=FALSE)
write.csv(cbind(1/2*(h.XXtw[[3]]$breaks[-length(h.XXtw[[3]]$breaks)]+h.XXtw[[3]]$breaks[-1]),h.XXtw[[3]]$density), file = "Xwt3pdf.csv", row.names=FALSE)
write.csv(cbind(1/2*(h.XXtw[[5]]$breaks[-length(h.XXtw[[5]]$breaks)]+h.XXtw[[5]]$breaks[-1]),h.XXtw[[5]]$density), file = "Xwt5pdf.csv", row.names=FALSE)


x11()
#dev.new() 
plot(1/2*(h.XX$breaks[-length(h.XX$breaks)]+h.XX$breaks[-1]),h.XX$density, main="Lifetime X pdf", xlab="x", ylab="f_X(x)", xlim=c(0,200000), ylim=c(0,3E-5), type="l", col="black")
#matplot(uxx,fTteor42,, pch=2, lty=1, add=T, type="l", col="blue")
#matplot(1/2*(TTh45$breaks[-length(TTh45$breaks)]+TTh45$breaks[-1]),TTh45$density,, pch=3, lty=1, add=T, col="red", type="p")
dev.print(device=pdf, "LifetimePdf.pdf")

x11()
#dev.new() 
plot(1/2*(h.XXtw[[1]]$breaks[-length(h.XXtw[[3]]$breaks)]+h.XXtw[[3]]$breaks[-1]),h.XXtw[[3]]$density, main="Pdf of RUL Xt|wt of liner 3", xlab="x", ylab="f_Xt|wt(x)", xlim=c(0,200000), ylim=c(0,3E-5), type="l", col="black")
#matplot(uxx,fTteor42,, pch=2, lty=1, add=T, type="l", col="blue")
#matplot(1/2*(TTh45$breaks[-length(TTh45$breaks)]+TTh45$breaks[-1]),TTh45$density,, pch=3, lty=1, add=T, col="red", type="p")
dev.print(device=pdf, "RUL3Pdf.pdf")

x11()
#dev.new() 
plot(1/2*(h.XXtw[[1]]$breaks[-length(h.XXtw[[5]]$breaks)]+h.XXtw[[5]]$breaks[-1]),h.XXtw[[5]]$density, main="Pdf of RUL Xt|wt of liner 5", xlab="x", ylab="f_Xt|wt(x)", xlim=c(0,200000), ylim=c(0,3E-5), type="l", col="black")
#matplot(uxx,fTteor42,, pch=2, lty=1, add=T, type="l", col="blue")
#matplot(1/2*(TTh45$breaks[-length(TTh45$breaks)]+TTh45$breaks[-1]),TTh45$density,, pch=3, lty=1, add=T, col="red", type="p")
dev.print(device=pdf, "RUL5Pdf.pdf")

# -----------------------
# SAVE Workspace
# -----------------------
save.image(file="BTGPanalysis.RData")




