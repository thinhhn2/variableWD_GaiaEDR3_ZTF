data.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/16 Pulsating White Dwarfs Candidates with Gaia Time-series Analysis/All Data - Numbered as in the Report/ZTF DR5"
plot.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/16 Pulsating White Dwarfs Candidates with Gaia Time-series Analysis/ZTF DR5 Plots/Plots with Guidry cleaning method"
name.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/16 Pulsating White Dwarfs Candidates with Gaia Time-series Analysis/ZTF DR5 Plots/ZTF g freq 16 candidates.csv"
name.data = read.csv(name.file,colClasses=c("name"="character"))
name = name.data$name

res = data.freq.gZTF(name,data.dir,plot.dir)

data.freq.gZTF = function(name,data.dir,plot.dir) { 
  library(FITSio)
  library(Rfast)
  frequency.data.1 = c()
  frequency.data.2 = c()
  frequency.data.3 = c()
  frequency.data.4 = c()
  frequency.data.5 = c()
  power.data.1 = c()
  power.data.2 = c()
  power.data.3 = c()
  power.data.4 = c()
  power.data.5 = c()
  sig.level.data = c()
  p.value.data = c()
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
    
    filename.lc = sprintf("%s Light Curve Full g band.png",n)
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
    
    timeG = timeG[catflag==0 & abs(sharp)<0.25 & G < limitmag - 1.0]
    err = err[catflag==0 & abs(sharp)<0.25 & G < limitmag - 1.0]
    G = G[catflag==0 & abs(sharp)<0.25 & G < limitmag - 1.0]
    
    filename.lc = sprintf("%s Light Curve g band.png",n)
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
    timeG.lsp = timeG - mean(timeG)
    G.lsp = G - mean(G)
    
    
    filename.lsp = sprintf("%s LSP g band.png",n)
    png(file=filename.lsp,width=700,height=520)
    par(mar=c(6,6,5,3)+.1)
    lsp.data = lsp(G.lsp,timeG.lsp,from=0,to=1000,ofac=8,alpha=0.01,plot=TRUE,xlab="",ylab="",main="",cex.lab=2.5,cex.axis=2)
    frequency = lsp.data$peak.at[1]
    title(sprintf("%s ZTF Periodogram in g Band \nPeak at %g [1/day]",n,frequency),cex.main=2.5)
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
    
    filename.pc = sprintf("%s Phase Curve g band.png",n)
    png(file=filename.pc,width=700,height=520)
    par(mar=c(6,6,5,3)+.1)
    xlab = xlab
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
    index.2 = nth(power.list,2,descending=TRUE,index.return=TRUE)
    max.power.2 = power.list[index.2]
    frequency.2 = freq.list[index.2]
    index.3 = nth(power.list,3,descending=TRUE,index.return=TRUE)
    max.power.3 = power.list[index.3]
    frequency.3 = freq.list[index.3]
    index.4 = nth(power.list,4,descending=TRUE,index.return=TRUE)
    max.power.4 = power.list[index.4]
    frequency.4 = freq.list[index.4]
    index.5 = nth(power.list,5,descending=TRUE,index.return=TRUE)
    max.power.5 = power.list[index.5]
    frequency.5 = freq.list[index.5]
    sig.level = lsp.data$sig.level #level for the power to be significant according to the lsp
    p.value = lsp.data$p.value #probability that the peak occurs by chance
    num.obs = lsp.data$n 
    
    frequency.data.1 = c(frequency.data.1,frequency.1)
    power.data.1 = c(power.data.1,max.power.1)
    frequency.data.2 = c(frequency.data.2,frequency.2)
    power.data.2 = c(power.data.2,max.power.2)
    frequency.data.3 = c(frequency.data.3,frequency.3)
    power.data.3 = c(power.data.3,max.power.3)
    frequency.data.4 = c(frequency.data.4,frequency.4)
    power.data.4 = c(power.data.4,max.power.4)
    frequency.data.5 = c(frequency.data.5,frequency.5)
    power.data.5 = c(power.data.5,max.power.5)
    sig.level.data = c(sig.level.data,sig.level)
    num.obs.data = c(num.obs.data,num.obs)
    p.value.data = c(p.value.data,p.value)
    time.max.data = c(time.max.data,max(timeG))
    time.min.data = c(time.min.data,min(timeG))
    name.data = c(name.data,n)
    band.data = c(band.data,"g")
    result.freq.power = data.frame("name"=name.data,"band"=band.data,"time.max"=time.max.data,"time.min"=time.min.data,"sig.level"=sig.level.data,"num.obs"=num.obs.data,"p.value"=p.value.data,"frequency.1"=frequency.data.1,"max.power.1"=power.data.1,"frequency.2"=frequency.data.2,"max.power.2"=power.data.2,"frequency.3"=frequency.data.3,"max.power.3"=power.data.3,"frequency.4"=frequency.data.4,"max.power.4"=power.data.4,"frequency.5"=frequency.data.5,"max.power.5"=power.data.5)
  }
  return(result.freq.power)
}
