file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/Reference for Periodic Percentage/Data/r band timeseries.csv"
file.query = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/Reference for Periodic Percentage/Data/r band query result.csv"
plot.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/Reference for Periodic Percentage/Plots"

ztf.data = read.csv(file=file,header=TRUE,stringsAsFactor=TRUE)
ztf.query = read.csv(file=file.query,header=TRUE,stringsAsFactor=TRUE) 
levels = levels(factor(ztf.data$oid))

setwd(plot.dir)

for (i in 1:length(levels)) {
  G = ztf.data$mag[ztf.data$oid == levels[i]]
  timeG = ztf.data$mjd[ztf.data$oid == levels[i]]
  err = ztf.data$magerr[ztf.data$oid == levels[i]]
  n = ztf.query$cntr_01[ztf.query$oid == levels[i]]
  
  filename.lc = sprintf("%d Light Curve Full r band.png",n)
  png(file=filename.lc,width=700,height=520)
  par(mar=c(6,6,5,3)+.1)
  plot(timeG,G,ylim=c(max(G)+err[G==max(G)],min(G)-err[G==min(G)]),xlab="",ylab="",col="red",cex=1.6,cex.lab=2.5,cex.axis=2,cex.main=2.5,pch=16,main= sprintf("%s ZTF Full Light Curve\n in r band",n))
  title(xlab="JD - 2400000.5",cex.lab=2.5,line=4)
  title(ylab="ZTF r mag",cex.lab=2.5,line=3.5)
  arrows(timeG, G-err, timeG, G+err, length=0.05, angle=90, code=3, col="red")
  box(lwd=2)
  axis(1,lwd=2,cex.axis=2)
  axis(2,lwd=2,cex.axis=2)
  dev.off()
  
  len = length(timeG)
  limit = len    #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
  for (i in 1:(len-1)){
    if(timeG[i+1]-timeG[i] >150) {
      limit = i
    }
  }
  timeG.left = timeG[1:limit]
  G.left = G[1:limit]
  err.left = err[1:limit]
  timeG.right = timeG[(limit+1):len]
  G.right = G[(limit+1):len]
  err.right = err[(limit+1):len]
  if (length(timeG.left)>length(timeG.right)){
    timeG = timeG.left
    G = G.left
    err = err.left }
  if (length(timeG.right)>length(timeG.left)){
    timeG = timeG.right
    G = G.right
    err = err.right }
  
  len = length(timeG)
  limit2 = len   #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
  for (i in 1:(len-1)){
    if(timeG[i+1]-timeG[i] >150) {
      limit2 = i  #make a new variable limit2 because limit already has a value
    }
  }
  timeG.left = timeG[1:limit2]
  G.left = G[1:limit2]
  err.left = err[1:limit2]
  timeG.right = timeG[(limit2+1):len]
  G.right = G[(limit2+1):len]
  err.right = err[(limit2+1):len]
  if (length(timeG.left)>length(timeG.right)){
    timeG = timeG.left
    G = G.left
    err = err.left }
  if (length(timeG.right)>length(timeG.left)){
    timeG = timeG.right
    G = G.right
    err = err.right }
  
  iqr = IQR(G)
  sd = iqr*20/27
  median = median(G)
  upper = median + 3*sd
  lower = median - 3*sd
  timeG = timeG[G<upper & G>lower]
  err = err[G<upper & G>lower]
  G = G[G<upper & G>lower]
  
  filename.lc = sprintf("%d Light Curve r band.png",n)
  png(file=filename.lc,width=700,height=520)
  par(mar=c(6,6,5,3)+.1)
  plot(timeG,G,ylim=c(max(G)+err[G==max(G)],min(G)-err[G==min(G)]),xlab="",ylab="",cex=1.6,cex.lab=2.5,cex.axis=2,cex.main=2.5,pch=16,col="red",main= sprintf("%s ZTF Light Curve\n in r band",n))
  title(xlab="JD - 2400000.5",cex.lab=2.5,line=4)
  title(ylab="ZTF r mag",cex.lab=2.5,line=3.5)
  arrows(timeG, G-err, timeG, G+err, length=0.05, angle=90, code=3,col="red")
  box(lwd=2)
  axis(1,lwd=2,cex.axis=2)
  axis(2,lwd=2,cex.axis=2)
  dev.off()
  
  library("bit64")
  library("lomb")
  timeG.lsp = timeG - mean(timeG)
  G.lsp = G - mean(G)
  
  
  filename.lsp = sprintf("%d LSP r band.png",n)
  png(file=filename.lsp,width=700,height=520)
  par(mar=c(6,6,5,3)+.1)
  lsp.data = lsp(G.lsp,timeG.lsp,from=0,to=500,ofac=8,alpha=0.01,plot=TRUE,xlab="",ylab="",main="",cex.lab=2.5,cex.axis=2)
  frequency = lsp.data$peak.at[1]
  title(sprintf("%s ZTF Periodogram in r Band \nPeak at %g [1/day]",n,frequency),cex.main=2.5)
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
  
  filename.pc = sprintf("%d Phase Curve in r band.png",n)
  png(file=filename.pc,width=700,height=520)
  par(mar=c(6,6,5,3)+.1)
  plot(phase,G,ylim=c(max(G.1)+err[G.1==max(G.1)],min(G.1)-err[G.1==min(G.1)]),cex.main=2.5,pch=16,col="red",main=sprintf("%s ZTF Folded Curve \nin r Band",n),ylab="",xlab="",cex=1.6,cex.lab=2.5,cex.axis=2)
  title(xlab=xlab,cex.lab=2.5,line=4)
  title(ylab="ZTF r mag",cex.lab=2.5,line=3.5)
  arrows(phase, G-err, phase, G+err, length=0.05, angle=90, code=3, col="red")
  box(lwd=2)
  axis(1,lwd=2,cex.axis=2)
  axis(2,lwd=2,cex.axis=2)
  
  dev.off()
}
  