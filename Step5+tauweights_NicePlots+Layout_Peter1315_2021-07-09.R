# Nice Plots for exit time
# SRC 2020-06-09, updated 2020-08-01

rm(list = ls())
graphics.off()

# Load results of Step 4

# Save results for further plotting
#save(Tstep,Xvar,D1s$x,D1s$y,D2,sigma,xeq,xvec,drift,diff,negPF,
#     xvec.ep,EPF,x,wts,ETL,ETR,xL,wtsL,xR,wtsR,
#     file='Results_tau_plots_Peter1315.Rdata')

load(file='Results_Step4ET_Peter1315.Rdata')

# Save data for Matlab
#stdlevel = as.data.frame(cbind(Tstep,Xvar))
#write.csv(stdlevel,file='Enrichment_Peter.csv')

# Calculate mean exit time for left basin ---------------------------------------------
nL = length(ETL[,1])
nR = length(ETR[,1])
nboth =nR + nL

#
meanETl = sum(ETL[1:nL,2]*wts[1:nL])/sum(wts[1:nL])/12
print('',quote=F)
print('Mean ET (hours) for left basin',quote=F)
print(meanETl)
print('-----------------------------------------------',quote=F)

# save axes
xL = x[1:nL]
wtsL = wts[1:nL]

# Calculate weighted average ET for right basin ===========================================================

#
meanETr = sum(ETR[1:nR,2]*wts[(nL):(nboth-1)])/sum(wts[(nL):(nboth-1)])/12
print('',quote=F)
print('Mean ET (hours) for right basin',quote=F)
print(meanETr)
print('-----------------------------------------------',quote=F)

xR = x[nL:(nboth-1)]
wtsR = wts[nL:(nboth-1)]
# decdoy is used for converting time to different units, if needed
decdoy = Tstep

# 5 panel version ----------------------------------------------------
windows(height=15,width=5)
par(mfrow=c(5,1),mar=c(3.9, 4.5, 0.5, 1) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
#
plot(decdoy,Xvar,type='l',col='black',lwd=1,xlab='Year',
     ylab='Chlorophyll')
#
par(mar=c(2, 4.5, 0, 1) + 0.1)
#
yrange=range(c(D1s$y,sigmas$y),na.rm=T)
xrange=range(D1s$x)
plot(D1s$x,D1s$y,ylim=yrange,type='l',col='black',lwd=2,xlab=' ', #xlab='Chlorophyll',
     ylab='Drift or Diffusion')
points(D1s$x,sigmas$y,type='l',col='black',lty=2,lwd=2)
abline(h=0,lwd=1,lty=2,col='gray')
legend(x=-1,y=0.28,legend=c('diff. (as s.d.)','drift'),lty=c(2,1),lwd=c(2,2),seg.len=4,
       col=c('black','black'),cex=1,text.font=2,bty='n')
#
plot(xvec.ep[2:100],EPF,type='l',lwd=2,col='black',xlim=xrange,xlab=' ', #xlab='Chlorophyll',
     ylab='Effective Potential')
#
xrange=range(c(ETL[,1],ETR[,1],D1s$x))
yrange=range(c(ETL[,2],ETR[,2]))/12
plot(ETL[,1],ETL[,2]/12,xlim=xrange,ylim=yrange,type='l',col='black',lwd=2,xlab=' ', #xlab='Chlorophyll',
     ylab='Exit Time, hours')
points(ETR[,1],ETR[,2]/12,type='l',col='black',lwd=2)
abline(v=xeq[2],lty=3,lwd=2)
#
par(mar=c(3.9, 4.5, 0, 1) + 0.1)
#
plot(x,wts,type='l',lwd=2,xlim=xrange,col='black',xlab='Chlorophyll',ylab='Density')

# Try Layout
m.plot = rbind(c(1,1),c(2,3),c(4,5))
print(m.plot)

windows(height=8,width=6)
layout(m.plot)
par(mar=c(3.9, 4.5, 0.5, 1) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
#
plot(decdoy,Xvar,type='l',col='black',lwd=1,xlab='Year',
     ylab='Chlorophyll Level')
text(x=2013.1,y=5.5,'A',cex=1.4,font=2)
#
par(mar=c(2, 4.5, 0, 1) + 0.1)
#
yrange=range(c(D1s$y,sigmas$y),na.rm=T)
xrange=range(D1s$x)
plot(D1s$x,D1s$y,ylim=yrange,type='l',col='black',lwd=2,xlab=' ', #xlab='Chlorophyll',
     ylab='Drift or Diffusion')
points(D1s$x,sigmas$y,type='l',col='black',lty=2,lwd=2)
abline(h=0,lwd=1,lty=2,col='gray')
legend(x=-7,y=0.3,legend=c('diffusion (as s.d.)','drift'),lty=c(2,1),lwd=c(2,2),seg.len=4,
       col=c('black','black'),cex=1,text.font=2,bty='n')
text(x=6,y=0.55,'B',cex=1.4,font=2)
#legend(x=-1,y=0.28,legend=c('diff. (as s.d.)','drift'),lty=c(2,1),lwd=c(2,2),seg.len=4,
#      col=c('black','black'),cex=1,text.font=2,bty='n')
#
plot(xvec.ep[2:100],EPF,type='l',lwd=2,col='black',xlim=xrange,xlab=' ', #xlab='Chlorophyll',
     ylab='Effective Potential')
text(x=6.6,y=-5.2,'C',cex=1.4,font=2)
#
par(mar=c(3.9, 4.5, 0, 1) + 0.1)
xrange=range(c(ETL[,1],ETR[,1],D1s$x))
yrange=range(c(ETL[,2],ETR[,2]))/12
plot(ETL[,1],ETL[,2]/12,xlim=xrange,ylim=yrange,type='l',col='black',lwd=2,xlab='Chlorophyll',
     ylab='Exit Time, hours')
points(ETR[,1],ETR[,2]/12,type='l',col='black',lwd=2)
abline(v=xeq[2],lty=3,lwd=2)
text(x=-4,y=5,paste('mean = ',round(meanETl,0)),cex=1.4,font=2)
text(x=4,y=5,paste('mean = ',round(meanETr,0)),cex=1.4,font=2)
text(x=5.9,y=22,'D',cex=1.4,font=2)
#
plot(x,wts,type='l',lwd=2,xlim=xrange,col='black',xlab='Chlorophyll',ylab='Density')
text(x=5.9,y=0.003,'E',cex=1.4,font=2)

# Try Layout + Color =================================================================

m.plot = rbind(c(1,1),c(2,3),c(4,5))
print(m.plot)

windows(height=8,width=6)
layout(m.plot)
par(mar=c(3.9, 4.5, 0.5, 1) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
#
plot(decdoy,Xvar,type='l',col='forestgreen',xaxt='n',lwd=1,xlab='Year',
     ylab='Chlorophyll')
axis(side = 1, pretty(Tstep, n = 4))
abline(h=xeq[2],col='darkred',lwd=2)
abline(v=c(2013:2016),lty=2,lwd=2,col='darkgray')
text(x=2013.1,y=7,'A',cex=1.5,font=2)

#
par(mar=c(2, 4.5, 0, 1) + 0.1)
#
yrange=range(c(D1s$y,sigmas$y),na.rm=T)
xrange=range(D1s$x)
plot(D1s$x,D1s$y,ylim=yrange,type='l',col='black',lwd=2,xlab=' ', 
     ylab='Drift or Diffusion')
points(D1s$x,sigmas$y,type='l',col='red',lty=2,lwd=2)
abline(h=0,lwd=1,lty=2,col='gray')
text(x=6.5,y=0.19,'B',cex=1.5,font=2)
legend(x=-7,y=0.1,legend=c('diffusion (as s.d.)','drift'),lty=c(2,1),lwd=c(2,2),seg.len=4,
       col=c('red','black'),cex=1,text.font=2,bty='n')
#
plot(xvec.ep[1:99],EPF,type='l',lwd=2,col='black',xlim=xrange,xlab=' ', 
     ylab='Effective Potential')
text(x=6,y=-8.6,'C',cex=1.5,font=2)
#
par(mar=c(3.9, 4.5, 0, 1) + 0.1)
xrange=range(c(ETL[,1],ETR[,1],D1s$x))
yrange=range(c(ETL[,2],ETR[,2]))/12
plot(ETL[,1],ETL[,2]/12,xlim=xrange,ylim=yrange,type='l',col='blue',lwd=2,xlab='Chlorophyll',
     ylab='Exit Time, hours')
points(ETR[,1],ETR[,2]/12,type='l',col='forestgreen',lwd=2)
abline(v=xeq[2],lty=3,lwd=2)
text(x=6.5,y=200,'D',cex=1.5,font=2)
text(x=-4.5,y=5,paste('mean = ',round(meanETl,0)),cex=1.4,font=2)
text(x=4.5,y=5,paste('mean = ',round(meanETr,0)),cex=1.4,font=2)
#
plot(x,wts,type='l',lwd=1,xlim=xrange,col='black',xlab='Chlorophyll',ylab='Density')
#polygon(c(xL,xL[length(xL)],xL[1]),c(wtsL,wtsL[1],wtsL[1]),col='skyblue',border=NA)
#polygon(c(xR,xR[1],xR[1]),c(wtsR,wtsR[length(wtsR)],wtsR[1]),col='lightgreen',border=NA)
polygon(x,wts,col='gray',border=NA)
text(x=6,y=0.01,'E',cex=1.5,font=2)
points(x,wts,type='l',lwd=2,xlim=xrange,col='black')

# Layout similar to Fig 1 of Science paper -------------------------------------------
#xrange = range(c(Xvar,xvec,D1s$x,ETL[,1],ETR[,1]),na.rm=T)
xrange = c(-7,7)
m.plot = cbind(c(1,1,2,3,4,5,6),c(1,1,1,1,1,1,1))
print(m.plot)
windows(width=6,height=10)
layout(m.plot)
par(mfrow=c(6,1),mar=c(1.2, 4.5, 1.1, 1.5) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
#
plot(Xvar,decdoy,type='l',lwd=1,col='green',xlim=xrange,xlab=' ',ylab='Year')
abline(v=xeq[2],col='black',lty=2,lwd=2)
text(x=-7,y=2015.8,'A',cex=1.5,font=2)
#
par(mar=c(1.2, 4.5, 0.6, 1.5) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
y.xeq = c(0,0,0)
pch.xeq = c(19,21,19)
plot(D1x,D1y,type='l',lwd=2,col='blue',xlim=xrange,ylim=c(-0.01,0.01),
     xlab=' ', ylab='Drift')
abline(h=0,lty=1,lwd=1,col='black')
points(c(xeq[1:2],xeq[3]-0.1),y.xeq,type='p',pch=pch.xeq,cex=2,lwd=2,col='black')
text(x=7,y=0.007,'B',cex=1.5,font=2)
#
plot(D2s$x,D2s$y,type='l',lwd=2,col='red',xlim=xrange,xlab=' ', ylab='Diffusion')
text(x=7,y=0.02,'C',cex=1.5,font=2)
#
# CHOOSE A POTENTIAL, COMMENT-OUT THE ALTERNATIVE
# Traditional (deterministic) potential
#plot(xvec,negPF,type='l',lwd=2,col='blue',xlim=xrange,xlab=' ',ylab='Potential')
#text(x=4.7,y=0,'D',cex=1.5,font=2)
#
# Effective (stochastic) potential
plot(xvec.ep[2:100],EPF,type='l',lwd=2,col='blue',xlim=xrange,xlab=' ',
     ylab='Effective Potential')
text(x=6,y=-8.8,'D',cex=1.5,font=2)
#
plot(x,wts,type='l',lwd=1,xlim=xrange,col='black',xlab=' ',ylab='Density')
#polygon(c(-7,xL,xL[length(xL)],xL[1]),c(0,wtsL,wtsL[1],wtsL[1]),col='skyblue',border=NA)
#polygon(c(0.5,xR,xR[1],xR[1]),c(0,wtsR,wtsR[length(wtsR)],wtsR[1]),col='lightgreen',border=NA)
polygon(x,
        wts,
        col='lightgreen',border=NA)
points(x,wts,type='l',lwd=2,xlim=xrange,col='black')
text(x=7,y=0.009,'E',cex=1.5,font=2)
#
par(mar=c(3.8, 4.5, 0.6, 1) + 0.1,cex.axis=1.1,font.axis=2,
    cex.lab=1.4,font.lab=2)
plot(ETL[,1],ETL[,2]/12,xlim=xrange,ylim=yrange,type='l',col='blue',lwd=2,
     xlab='Phycocyanin level',ylab='Exit Time, h')
points(ETR[,1],ETR[,2]/12,type='l',col='forestgreen',lwd=2)
abline(v=xeq[2],col='black',lty=2,lwd=2)
text(x=7,y=199,'F',cex=1.5,font=2)
text(x=-4.5,y=20,paste('mean = ',round(meanETl,0)),cex=1.4,font=2)
text(x=4.5,y=20,paste('mean = ',round(meanETr,0)),cex=1.4,font=2)
