# Test DDbintau
# SRC 2021-02-16

rm(list = ls())
graphics.off()

library(stats)
library(moments)

source('DDbintau+D4.r')

# save results
# Save post-processed DLM data
# X.raw is original series, X.z is z score, X.rawmean and X.rawsd are mean & sd for z score
# nl is DLM lags, delta is DLM delta, X.dlm is input to DLM, Tstep is time,
# local equilibria are X.eq and SD.eq, with ratio z.eq
# level (intercept) measures are level or stdlevel (for level/sd)
# Yyhat, B.ests, B.sd, errvar are DLM outputs
#save(nl,delta,X.raw,X.z,X.rawmean,X.rawsd,X.dlm,Tstep,X.eq,SD.eq,level,stdlevel,Z.eq,
#     Yyhat,B.ests,B.sd,errvar,file=Fname)

# Choose dataset 
# Peter 1315
#load(file='Peter1315_DLMresult.Rdata')
#fname=c('DDbintau_Peter1315.Rdata')
#plot.title=c('Peter, Enrich')
#
# Tuesday 1315
load(file='Tuesday1315_DLMresult.Rdata')
fname=c('DDbintau_Tuesday1315.Rdata')
plot.title=c('Tuesday, Enrich')

# Histograms of DLM results
s.level = level/stdlevel # s(b,t) = b(t)/(b(t)/s(b,t))
windows(height=12,width=6)
par(mfrow=c(4,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
hist(X.z,breaks=200,col='skyblue')
hist(level,breaks=200,col='skyblue')
plot(level,s.level,log='y',type='p',pch=19,cex=0.3,xlab='b(t)',ylab='s(b,t)')
abline(h=1,lty=2)
hist(stdlevel,breaks=200,col='skyblue')

# Choose variate for Langevin model from DLM output
Xvar = na.omit(stdlevel) #[10:length(stdlevel)]
nxvar = length(Xvar)
Tstep = Tstep[1:length(Xvar)]  

# extract 2015 only
#T15 = which(Tstep >= 2015)
#nxvar = length(T15)
#Xvar = Xvar[T15]
#Tstep = Tstep[T15]
#fname=c('DDbintau_level_Peter2015.Rdata')

windows(width=10,height=4)
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep,Xvar,type='l',lwd=2,xlab='day',
     ylab='Xvar')
grid()

windows(width=10,height=4)
par(mfrow=c(1,1),mar=c(4, 4.2, 2, 2) + 0.1,cex.axis=1.5,cex.lab=1.5)
plot(Tstep[1:(nxvar-1)],diff(Xvar),type='l',lwd=2,xlab='day',
     ylab='diff(Xvar)')
grid()

# Set up for binning method
# DDbins = function(Xvar,bw,ntau,nbin)
title = c('Pigment variate')
ntau = 15 # include lag with minimum X2 from Chapman-Kolmogorov; see Matlab plot
nbin = 200 # usually 100 to 200
print(c('N per bin = ',nxvar/nbin),quote=F)
#bw <- 0.3*sd(Xvar)  # mesh method: try between 0.1*sd and 0.5*sd 
bw = 0.1*diff(range(Xvar)) # tie bw to range of time series for bins
# run function
DDout = DDbins(Xvar,bw,ntau,nbin)
# extract smoothed output
#outlist = list(D1s,D2s,sigmas,D4s,bin.mid,D1,D2,sigma,D4,SD) 
D1s = DDout[[1]]
D2s = DDout[[2]]
sigmas = DDout[[3]]
D4s = DDout[[4]]
bin.mid= DDout[[5]]
D1 = DDout[[6]]
D2 = DDout[[7]]
D4 = DDout[[9]]
SD = DDout[[10]]
SDall = DDout[[11]]

# look at results of DDbintau
windows()
par(mfrow=c(3,1),mar=c(4,4.5,2,2)+0.1,cex.axis=2,cex.lab=2)
plot(D1s$x,D1s$y,type='l',lwd=2,col='blue',xlab='state',ylab='D1')
abline(h=0,lty=2)
plot(sigmas$x,sigmas$y,type='l',lwd=2,col='red',xlab='state',ylab='sigma')
plot(D2s$x,D2s$y,type='l',lwd=2,col='red',xlab='state',ylab='D2')

# Find equilibria
sdrift = sign(D1s$y)
dsdrift = c(0,-diff(sdrift))
xeq = D1s$x[which(!dsdrift == 0)]
ixeq = which(!dsdrift == 0)  # indices of the equilibria

print('equilibria and their indices',quote=F)
print(xeq,quote=F)
print(ixeq,quote=F)

# Test one-step predictions

# Function for negative log likelihood
NLLfun = function(delta) {  # delta is residual 
  xN = length(delta)
  xmu = mean(delta,na.rm=T) # not needed for delta = x-xhat
  xsig = sd(delta,na.rm=T)
  # from Hilborn & Mangel eq 7.11 p. 137
  T1 = xN*(log(xsig) + 0.5*log(2*pi))
  T3 = (1/(2*xsig^2))*sum(delta^2)
  NLL = T1 + T3
  return(NLL)
}

# Write D1 and D2 as functions 
D1nona = na.omit(as.data.frame(D1s))
D1spline = smooth.spline(x=D1nona$x,y=D1nona$y)
D1fct = function(x) {
  #yhat = approx(x=avec, y=Drift.vec,xout=x,method='linear',rule=2)$y
  yhat=predict(D1spline,x)$y
  return(yhat)
}
D2nona = na.omit(as.data.frame(D2s))
D2spline = smooth.spline(x=D2nona$x,y=D2nona$y)
D2fct = function(x)  {
  #yhat = approx(x=avec, y=D2.from.Sig,xout=x,method='linear',rule=2)$y
  yhat=predict(D2spline,x)$y
  return(yhat)
}

# Make a model with constant D1 and use it in Langevin
D1flat = mean(D1nona$y,na.rm=T)
D2flat = mean(D2nona$y,na.rm=T)

Xvar0 = Xvar[1:(nxvar-1)]
Xvar1 = Xvar[2:nxvar]
NXvar0 = length(Xvar0)

# If needed, remove Xvar starting points that are out of bounds for D1 and D2
#x.range = range(c(D1nona$x,D2nona$x),na.rm=T)
#xmat.all = as.data.frame( cbind(Xvar0,Xvar1) )
#xmat.nona = na.omit(xmat.all)
#xmat.ok = subset(xmat.nona,subset=(xmat.nona > x.range[1] & xmat.nona < x.range[2] ))

#x0 = as.vector(xmat.ok$Xvar0)
#x1 = as.vector(xmat.ok$Xvar1)
#print(c('dim of data not NA and within range ',dim(xmat.ok)),quote=F)
#nxsim = length(x0)

# Constants used in Langevin 
dt = 1/288 # in days; 288 samples/day
dtnoise = sqrt(dt)
dW = rnorm(NXvar0,mean=0,sd=1)

# Run Langevin with flat D1
dev0 = as.vector(rep(0,NXvar0))
for(i in 2:NXvar0) {
  xx = Xvar0[i]
  sig = sqrt(2*D2fct(xx))
  sig.flat = sqrt(2*D2flat)
  #sig = ifelse(is.na(sig),sig.flat,sig)
  #dxhat = D1fct(xx)*dt + sig*dtnoise*dW[i] # observed D1
  dxhat = D1flat*dt + sig.flat*dtnoise*dW[i] # flat D1 & sigma
  #dxhat = D1flat*dt + sig*dtnoise*dW[i] # flat D1 & observed sigma
  x1hat = xx + dxhat
  dev0[i] = (Xvar1[i] - x1hat)
}

dev0nona = na.omit(dev0)
windows()
plot(dev0nona,type='l',xlab='Time',ylab='One-step error')

# variance
v.F = var(dev0nona)
print('Langevin one-step variance - flat D1',quote=F)
print(v.F)
# NLL
NLL.F = NLLfun(dev0nona)
print('Negative log likelihood & aic for Flat D1',quote=F)
print(NLL.F)
aic.F = 4 + 2*NLL.F
print(aic.F)

# variance of one-step predictions from fitted Langevin

dt = 1/288 # in days; 288 samples/day
dev = rep(0,NXvar0)
for(i in 2:NXvar0) {
  xx = Xvar0[i]
  sig = sqrt(2*D2fct(xx))
  dxhat = D1fct(xx)*dt + sig*dtnoise*dW[i]
  x1hat = xx + dxhat
  dev[i] = (Xvar1[i] - x1hat)
}

devnona = na.omit(dev)
windows()
plot(devnona,type='l',xlab='Time',ylab='One-step error')

# variance
v.L = var(devnona)
print('Langevin one-step variance, actual D1',quote=F)
print(v.L)
# NLL
NLL.L = NLLfun(devnona)
print('Negative log likelihood & aic for actual D1',quote=F)
print(NLL.L)
#aic.L = 4 + 2*NLL.L
#print(aic.L)

# delta AIC
print('delta NLL, observed - flat D1',quote=F)
print(NLL.L - NLL.F)

# save Langevin fit for next steps
save(Tstep,X.raw,Xvar,bin.mid,D1s,D2s,sigmas,xeq,file=fname)
