library("R2OpenBUGS")
library(coda)
require(devEMF)

getwd() -> workdir
suffnome <- "wlimGamma_bGamma_betaUnifVague_aUnifVague" # priors U (wlim here) and b 3-params Gamma, a and beta vague Unif

codpkg = FALSE # set to FALSE
DICval = FALSE

thin = 200 
burnin = 1.0*10^5
iteraz = 1.0*10^5 + burnin

nochains = 1 ### no. of chains [1-3]  ####

init1 <- list(b=1.0, beta=10.0, a=100, wlim=4.7) # initialization can be changed
init2 <- list(b=1.1, beta=9.0, a=150, wlim=4.6) # 
init3 <- list(b=1.2, beta=8.0, a=200, wlim=4.6) # 
if (nochains==3) {
  initial <-list(init1,init2,init3) # if 3 chain
  } else if (nochains==2) {
    initial <-list(init1,init2) # if 2 chain
} else initial <- list(init1) # if 1 chain

simnum = iteraz - burnin

Sys.time()->tstart
bayes.sim3 <- bugs("LinersData.txt", inits = initial, model.file = paste0("BoundedTransfGamma_",suffnome,".txt"), parameters = c("beta", "wlim", "a", "b"), n.chains = nochains, n.iter = iteraz, n.burnin = burnin, n.thin = thin, DIC = DICval, codaPkg = codpkg, bugs.seed=2)

Tbugs <- Sys.time()-tstart
cat(sprintf("\nTempo bayes.sim3 con %d catena/e: %f s\n\n",nochains,as.numeric(Tbugs, units = "secs")))

# Create a new dir for results, if it does not exist
if (file.exists(paste0("Results_",suffnome))){
  setwd(paste0("Results_",suffnome))
} else {
  dir.create(paste0("Results_",suffnome))
  setwd(paste0("Results_",suffnome))
}
  
# print(bayes.sim)  # quando codaPkg = FALSE
# plot(bayes.sim)   # quando codaPkg = FALSE


Sys.time()->tstart
if (codpkg==FALSE) {
if (nochains==3) {
  codaobject1 <- as.mcmc(bayes.sim3$sims.array[,1,]) # quando codaPkg = FALSE
  codaobject2 <- as.mcmc(bayes.sim3$sims.array[,2,]) # quando codaPkg = FALSE
  codaobject3 <- as.mcmc(bayes.sim3$sims.array[,3,]) # quando codaPkg = FALSE
} else if (nochains==2) {
  codaobject1 <- as.mcmc(bayes.sim3$sims.array[,1,]) # quando codaPkg = FALSE
  codaobject2 <- as.mcmc(bayes.sim3$sims.array[,2,]) # quando codaPkg = FALSE
} else codaobject1 <- as.mcmc(bayes.sim3$sims.array[,1,]) # quando codaPkg = FALSE
} else if (nochains==3) {  # quando codaPkg = TRUE
  codaobject1 <- read.bugs(bayes.sim3[[1]])
  codaobject2 <- read.bugs(bayes.sim3[[2]])
  codaobject3 <- read.bugs(bayes.sim3[[3]])
} else if (nochains==2) { # quando codaPkg = TRUE
  codaobject1 <- read.bugs(bayes.sim3[[1]])
  codaobject2 <- read.bugs(bayes.sim3[[2]])
} else codaobject1 <- read.bugs(bayes.sim3) # quando codaPkg = TRUE e una catena

#save.image(file="workspaceRunCODA.RData")

if (nochains==3) {
  mc.out1 <- codaobject1[1:simnum,]
  mc.out2 <- codaobject2[1:simnum,]
  mc.out3 <- codaobject3[1:simnum,]
} else if (nochains==2) {
  mc.out1 <- codaobject1[1:simnum,]
  mc.out2 <- codaobject2[1:simnum,]
} else mc.out1 <- codaobject1[1:simnum,]



print(" ")
print("Tempo CODA e mc.out1")
print(Sys.time()-tstart)

# VISUAL INSPECTION DIAGNOSTICS
print(" ")
print(" ")
print(" VISUAL inspection diagnostics")
print(" ")
print(" ")

# Sys.time()->tstart;
# print(" Summary of CODA object")
# print(" ")
if (codpkg==TRUE) {
  if (nochains==3) {
    print(summary(codaobject1))
    print(summary(codaobject2))
    print(summary(codaobject3))
  }
  if (nochains==2) {
    print(summary(codaobject1))
    print(summary(codaobject2))
  }
  if (nochains==1) print(summary(codaobject1))
}

print(" Summary of MC no.1")
print(" ")
print(bayes.sim3)

print(" ")
print("Tempo summary")
print(Sys.time()-tstart)

x11()

  if (nochains==3) {
      plot(codaobject1,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity1_",suffnome,".pdf"))
      x11()
      plot(codaobject2,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity2_",suffnome,".pdf"))
      x11()
      plot(codaobject3,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity3_",suffnome,".pdf"))
    }
 if (nochains==2) {
      plot(codaobject1,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity1_",suffnome,".pdf"))
      x11()
      plot(codaobject2,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity2_",suffnome,".pdf"))
    }
if (nochains==1) {
      plot(codaobject1,  auto.layout = TRUE)
      dev.print(device=pdf, paste0("TraceDensity1_",suffnome,".pdf"))
    }

# RUNNING MEAN PLOTS
print(" ")
print(" ")
print(" RUNNING MEAN PLOTS")
print(" ")
print(" ")

Sys.time()->tstart
if(codpkg==FALSE)
{
if (nochains==3) {
  nc=ncol(codaobject1)
  runmean <- matrix(NA, nrow=nrow(codaobject1), nc)
  x11() 
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject1[,ii])/seq(along=codaobject1[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean1_",suffnome,".pdf"))
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart)
  par(mfrow=c(1,1))
  nc=ncol(codaobject2)
  runmean <- matrix(NA, nrow=nrow(codaobject2), nc)
  x11()
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject2[,ii])/seq(along=codaobject2[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean2_",suffnome,".pdf"))
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart);
  par(mfrow=c(1,1))
  nc=ncol(codaobject3)
  runmean <- matrix(NA, nrow=nrow(codaobject3), nc)
   x11()
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject3[,ii])/seq(along=codaobject3[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean3_",suffnome,".pdf"))
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart);
  par(mfrow=c(1,1))
  
} else if (nochains==2) {
  nc=ncol(codaobject1)
  runmean <- matrix(NA, nrow=nrow(codaobject1), nc)
  x11()
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject1[,ii])/seq(along=codaobject1[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean1_",suffnome,".pdf"))
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart)
  par(mfrow=c(1,1))
  nc=ncol(codaobject2)
  runmean <- matrix(NA, nrow=nrow(codaobject2), nc)
  x11() 
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject2[,ii])/seq(along=codaobject2[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean2_",suffnome,".pdf"))
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart);
  par(mfrow=c(1,1))
} else {
  if (codpkg==FALSE)
    {
    nc=ncol(codaobject1)
  runmean <- matrix(NA, nrow=nrow(codaobject1), nc)
   x11()
  par(mfrow=c(2,2))
  for (ii in 1:nc){
    runmean[,ii] <- cumsum(codaobject1[,ii])/seq(along=codaobject1[,ii])
    plot(runmean[,ii])
    dev.print(device=pdf, paste0("RunMean1_",suffnome,".pdf"))
  } 
  } else {
    nc=ncol(codaobject1[[1]])
    runmean <- matrix(NA, nrow=nrow(codaobject1[[1]]), nc)
    x11() 
    par(mfrow=c(2,2))
    for (ii in 1:nc){
      runmean[,ii] <- cumsum(codaobject1[[1]][,ii])/seq(along=codaobject1[[1]][,ii])
      plot(runmean[,ii])
      dev.print(device=pdf, paste0("RunMean1_",suffnome,".pdf"))
  }
  
  }
  print(" ")
  print("Tempo running mean plots")
  print(Sys.time()-tstart)
  par(mfrow=c(1,1))  
}
}

# Autocorrelations plots
print(" ")
print(" ")
print(" Autocorrelations functions")
print(" ")

if (nochains==3) {
  print(autocorr(codaobject1, lags = c(0, 1, 2, 3, 4, 5, 6)))
  print(autocorr(codaobject2, lags = c(0, 1, 2, 3, 4, 5, 6)))
  print(autocorr(codaobject3, lags = c(0, 1, 2, 3, 4, 5, 6)))
} else if (nochains==2) {
  print(autocorr(codaobject1, lags = c(0, 1, 2, 3, 4, 5, 6)))
  print(autocorr(codaobject2, lags = c(0, 1, 2, 3, 4, 5, 6)))
} else print(autocorr(codaobject1, lags = c(0, 1, 2, 3, 4, 5, 6)))

print(" ")
print(" Autocorrelations plots")
print(" ")
print(" ")

Sys.time()->tstart

if (nochains==3) {
  x11()
  autocorr.plot(codaobject1, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr1_",suffnome,".pdf"))
  print(autocorr.diag(codaobject1))
  x11()
  autocorr.plot(codaobject2, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr2_",suffnome,".pdf"))
  print(autocorr.diag(codaobject2))
  x11()
  autocorr.plot(codaobject3, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr3_",suffnome,".pdf"))
  print(autocorr.diag(codaobject3))
} else if (nochains==2) {
  x11()
  autocorr.plot(codaobject1, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr1_",suffnome,".pdf"))
  print(autocorr.diag(codaobject1))
  x11()
  autocorr.plot(codaobject2, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr2_",suffnome,".pdf"))
  print(autocorr.diag(codaobject2))
} else {
  x11()
  autocorr.plot(codaobject1, lag.max=15)
  dev.print(device=pdf, paste0("AutoCorr1_",suffnome,".pdf"))
  print(autocorr.diag(codaobject1))
}

print(" ")
print("Tempo autocorr. plots")
print(Sys.time()-tstart)

# --------
# BROOKS-GELMAN-RUBIN diagnostics
# --------
if (nochains>1) {
  cat("\n\n GELMAN diagnostics\n\n")
  #mh.list <- mcmc.list(list(mc.out1, mc.out2, mc.out3))
  if(Sys.info()["sysname"] == "Windows") x11() else quartz()
  if(codpkg==TRUE){
    x11()
    gelman.plot(read.bugs(bayes.sim3))
  }
  else if (nochains==2){
    x11()
    gelman.plot(mcmc.list(list(codaobject1, codaobject2)))
  }
  else {
    x11()
    gelman.plot(mcmc.list(list(codaobject1, codaobject2, codaobject3)))
  }
  
  dev.print(device=pdf, paste0("GelmanPlot_",suffnome,".pdf"))
  
  if(codpkg==TRUE) print(gelman.diag(read.bugs(bayes.sim3)))
  else print(gelman.diag(bayes.sim3))
}

print(" ")
print(" ")
print(" GEWEKE diagnostics")
print(" ")
print(" ")

print(" ")
print(" ")
print(" Geweke diagnostics of MC no.1")
print(" ")

Sys.time()->tstart

if (nochains==3) {
  print(geweke.diag(codaobject1))
  print(geweke.diag(codaobject2))
  print(geweke.diag(codaobject3))
} else if (nochains==2) {
  print(geweke.diag(codaobject1))
  print(geweke.diag(codaobject2))
} else print(geweke.diag(codaobject1))


print(" ")
print("Tempo Geweke diagnostics")
print(Sys.time()-tstart)

x11()
if(codpkg==FALSE) {
  dev.new()
  dev.new()
  geweke.plot(bayes.sim3)
  } else {
  if(nochains==1) {
    dev.new()
    dev.new()
    geweke.plot(codaobject1)
    }
  else {
    dev.new()
    dev.new()
    geweke.plot(read.bugs(bayes.sim3))
  }
}

#geweke.plot(bayes.sim3)
dev.print(device=pdf, paste0("GewekePlot_",suffnome,".pdf"))



print(" ")
print(" ")
print(" RAFTERY-LEWIS diagnostics")
print(" ")
print(" ")


print(" ")
print(" ")
print(" Raftery-Lewis diagnostics of MC no.1")
print(" ")

Sys.time()->tstart

if (nochains==3) {
  print(raftery.diag(codaobject1, q = 0.025, r = 0.005, s = 0.95))
  print(raftery.diag(codaobject2, q = 0.025, r = 0.005, s = 0.95))
  print(raftery.diag(codaobject3, q = 0.025, r = 0.005, s = 0.95))
} else if (nochains==2) {
  print(raftery.diag(codaobject1, q = 0.025, r = 0.005, s = 0.95))
  print(raftery.diag(codaobject2, q = 0.025, r = 0.005, s = 0.95))
} else print(raftery.diag(codaobject1, q = 0.025, r = 0.005, s = 0.95))



print(" ")
print("Tempo Raftery-Lewis diagnostics")
print(Sys.time()-tstart)



print(" ")
print(" ")
print(" Heidelberger-Welch diagnostics")
print(" ")
print(" ")


print(" ")
print(" ")
print(" Heidelberger-Welch diagnostics of MC no.1")
print(" ")

Sys.time()->tstart

if (nochains==3) {
  print(heidel.diag(codaobject1))
  print(heidel.diag(codaobject2))
  print(heidel.diag(codaobject3))
} else if (nochains==2) {
  print(heidel.diag(codaobject1))
  print(heidel.diag(codaobject2))
} else print(heidel.diag(codaobject1))

print(" ")
print("Tempo Heidelberger-Welch diagnostics")
print(Sys.time()-tstart)


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

x11()
#dev.new()
h<-hist(mc.sel[,"b"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h$breaks[-length(h$breaks)]+h$breaks[-1]),h$density, main="Posterior pdf of b", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("bpostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h2<-hist(mc.sel[,"wlim"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h2$breaks[-length(h2$breaks)]+h2$breaks[-1]),h2$density, main="Posterior pdf of wlim", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("wlimpostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h3<-hist(mc.sel[,"a"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h3$breaks[-length(h3$breaks)]+h3$breaks[-1]),h3$density, main="Posterior pdf of a", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("apostpdf_",suffnome,".pdf"))

x11()
#dev.new()
h4<-hist(mc.sel[,"beta"],breaks = floor(simnum^.4), plot=F)
plot(1/2*(h4$breaks[-length(h4$breaks)]+h4$breaks[-1]),h4$density, main="Posterior pdf of beta", xlab="x", ylab="p(x)",type="b")
dev.print(device=pdf, paste0("betapostpdf_",suffnome,".pdf"))

Sys.time()->tstart;

cat("\n\nMean of a",mean(mc.sel[,'a']),'\n')
cat("\n\nMean of b",mean(mc.sel[,'b']),'\n')
cat("\n\nMean of wlim",mean(mc.sel[,'wlim']),'\n')
cat("\n\nMean of beta",mean(mc.sel[,'beta']),"\n\n\n")

print(" ")
print("0.90 equal-tails credibility interval for b")
print(" ")
qb <- quantile(mc.sel[,'b'],  probs = c(5, 95)/100)
print(qb)
print(" ")

print(" ")
print("0.90 equal-tails credibility interval for wlim")
print(" ")
qw <- quantile(mc.sel[,'wlim'],  probs = c(5, 95)/100)
print(qw)
print(" ")

print(" ")
print("0.90 equal-tails credibility interval for a")
print(" ")
qa <- quantile(mc.sel[,'a'],  probs = c(5, 95)/100)
print(qa)
print(" ")

print(" ")
print("0.90 equal-tails credibility interval for beta")
print(" ")
qbet <- quantile(mc.sel[,'beta'],  probs = c(5, 95)/100)
print(qbet)
print(" ")

print(" ")
print("Tempo calcolo medie e quantili")
print(Sys.time()-tstart);

## END Posterior and credible intervals


# HPD intervals
Sys.time()->tstart;

print(" ")
print("0.90 HPD interval for a")
print(" ")
Hpda <- HPDinterval(as.mcmc(mc.sel[,'a']),  prob = 0.9)
print(Hpda[1,])
print(" ")

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
print("0.90 HPD interval for beta")
print(" ")
Hpdbet <- HPDinterval(as.mcmc(mc.sel[,'beta']),  prob = 0.9)
print(Hpdbet[1,])
print(" ")


print(" ")
print("Tempo calcolo intervalli HPD")
print(Sys.time()-tstart);
print(" ")

#---------------
# SAVE Workspace
#---------------
save.image(file=paste0("BugsBTGPsamples_",suffnome,"_noDIC.RData"))
setwd(workdir)

