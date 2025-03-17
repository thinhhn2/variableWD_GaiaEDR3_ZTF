data.freq.gZTF = function(name,data.dir,plot.dir) { 
library(FITSio)
#library(Rfast)
frequency.data.1 = c()
power.data.1 = c()
num.obs.data = c()
time.max.data = c()
time.min.data = c()
name.data = c()
band.data = c()

for (n in name){
file = sprintf("%sg.fits",n)

setwd(data.dir)

Data = readFITS(file)
G = Data$col[[6]]
timeG = Data$col[[5]]
err = Data$col[[7]]
catflag = Data$col[[8]]
sharp = Data$col[[13]]
limitmag = Data$col[[18]]

setwd(plot.dir)

filename.lc = sprintf("%s Light Curve Full.png",n)
png(file=filename.lc,width=700,height=520)
par(mar=c(6,6,5,3)+.1)
plot(timeG,G,ylim=c(max(G)+err[G==max(G)],min(G)-err[G==min(G)]),xlab="",ylab="",cex=1.6,cex.lab=2.5,cex.axis=2,cex.main=2.5,pch=16,main= sprintf("%s ZTF Full Light Curve\n in g band",n))
title(xlab="JD - 2400000.5",cex.lab=2.5,line=4)
title(ylab="ZTF g mag",cex.lab=2.5,line=3.5)
arrows(timeG, G-err, timeG, G+err, length=0.05, angle=90, code=3)
box(lwd=2)
axis(1,lwd=2,cex.axis=2)
axis(2,lwd=2,cex.axis=2)
dev.off()

timeG = timeG[catflag == 0 & abs(sharp) < 0.25 & G < limitmag - 1.0]
err = err[catflag == 0 & abs(sharp) < 0.25 & G < limitmag - 1.0]
G = G[catflag == 0 & abs(sharp) < 0.25 & G < limitmag - 1.0]

filename.lc = sprintf("%s Light Curve.png",n)
png(file=filename.lc,width=700,height=520)
par(mar=c(6,6,5,3)+.1)
plot(timeG,G,ylim=c(max(G)+err[G==max(G)],min(G)-err[G==min(G)]),xlab="",ylab="",cex=1.6,cex.lab=2.5,cex.axis=2,cex.main=2.5,pch=16,main= sprintf("%s ZTF Light Curve\n in g band",n))
title(xlab="JD - 2400000.5",cex.lab=2.5,line=4)
title(ylab="ZTF g mag",cex.lab=2.5,line=3.5)
arrows(timeG, G-err, timeG, G+err, length=0.05, angle=90, code=3)
box(lwd=2)
axis(1,lwd=2,cex.axis=2)
axis(2,lwd=2,cex.axis=2)
dev.off()

library("bit64")
library("lomb")
#timeG.lsp = timeG - mean(timeG)
timeG.lsp = timeG
G.lsp = G - mean(G)


filename.lsp = sprintf("%s LSP.png",n)
png(file=filename.lsp,width=700,height=520)
par(mar=c(6,6,5,3)+.1)
lsp.data = lsp(G.lsp,timeG.lsp,from=0,to=400,ofac=8,alpha=0.01,plot=FALSE)
frequency = lsp.data$peak.at[1]
plot(lsp.data$scanned, lsp.data$power, type="l",xlab="",ylab="",main="",cex.lab=2.5,cex.axis=2)
title(sprintf("%s Periodogram g \nPeak at %g [1/day]",n,frequency),cex.main=2.5)
title(xlab="frequency[1/day]",cex.lab=2.5,line=4)
title(ylab="normalised power",cex.lab=2.5,line=3.5)
box(lwd=2)
axis(1,lwd=2,cex.axis=2)
axis(2,lwd=2,cex.axis=2)
dev.off()

period = 1/frequency
period.sec.G = period*24*60*60
phase1 = (timeG%%period)/period
phase2 = phase1 + 1
phase = c(phase1,phase2)
G.1 = G
G = c(G,G)

if(period.sec.G<60){xlab = sprintf("phase (Period = %.3f secs)",period.sec.G)}
if(period.sec.G>=60 & period.sec.G<3600){xlab = sprintf("phase (Period = %.3f mins)",period.sec.G/60)}
if(period.sec.G>=3600 & period.sec.G<86400){xlab = sprintf("phase (Period = %.3f hours)",period.sec.G/3600)}
if(period.sec.G>=86400 & period.sec.G<31536000){xlab = sprintf("phase (Period = %.3f days)",period.sec.G/86400)}
if(period.sec.G>=31536000){xlab = sprintf("phase (Period = %.3f years)",period.sec.G/31536000)}

filename.pc = sprintf("%s Phase Curve.png",n)
png(file=filename.pc,width=700,height=520)
par(mar=c(6,6,5,3)+.1)
plot(phase,G,ylim=c(max(G.1)+err[G.1==max(G.1)],min(G.1)-err[G.1==min(G.1)]),cex.main=2.5,pch=16,main=sprintf("%s ZTF Folded Curve \nin g Band",n),ylab="",xlab="",cex=1.6,cex.lab=2.5,cex.axis=2)
title(xlab=xlab,cex.lab=2.5,line=4)
title(ylab="ZTF g mag",cex.lab=2.5,line=3.5)
arrows(phase, G-err, phase, G+err, length=0.05, angle=90, code=3)
box(lwd=2)
axis(1,lwd=2,cex.axis=2)
axis(2,lwd=2,cex.axis=2)

dev.off()

freq.list = lsp.data$scanned
power.list = lsp.data$power

max.power.1 = lsp.data$peak
frequency.1 = lsp.data$peak.at[1]
num.obs = lsp.data$n 

frequency.data.1 = c(frequency.data.1,frequency.1)
power.data.1 = c(power.data.1,max.power.1)
num.obs.data = c(num.obs.data,num.obs)
time.max.data = c(time.max.data,max(timeG))
time.min.data = c(time.min.data,min(timeG))
name.data = c(name.data,n)
band.data = c(band.data,"g")
result.freq.power = data.frame("name"=name.data,"band"=band.data,"time.max"=time.max.data,"time.min"=time.min.data,"num.obs"=num.obs.data,"frequency.1"=frequency.data.1,"max.power.1"=power.data.1)
}
return(result.freq.power)
}


