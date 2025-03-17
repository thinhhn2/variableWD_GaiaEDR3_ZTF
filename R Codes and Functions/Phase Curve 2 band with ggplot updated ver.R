library(FITSio)
library(ggplot2)
library(cowplot)

data.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/DQV/Data"
plot.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/DQV/Plots/Periodic 2-band phase curves"

setwd(data.dir)

sourceid = "1191504471436192512"
frequency = 7.68743
alpha = 0.5 #factor for setting the y-axis limit

file.g = paste(sourceid,"g.fits",sep="")
file.r = paste(sourceid,"r.fits",sep="")

Data.g = readFITS(file.g)
Data.r = readFITS(file.r)

G = Data.g$col[[6]]
timeG = Data.g$col[[5]]
errG = Data.g$col[[7]]
R = Data.r$col[[6]]
timeR = Data.r$col[[5]]
errR = Data.r$col[[7]]

#Cleaning algorithm
len = length(timeG)
limit = len    #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
for (i in 1:(len-1)){
  if(timeG[i+1]-timeG[i] >150) {
    limit = i
  }
}
timeG.left = timeG[1:limit]
G.left = G[1:limit]
errG.left = errG[1:limit]
timeG.right = timeG[(limit+1):len]
G.right = G[(limit+1):len]
errG.right = errG[(limit+1):len]
if (length(timeG.left)>length(timeG.right)){
  timeG = timeG.left
  G = G.left
  errG = errG.left }
if (length(timeG.right)>length(timeG.left)){
  timeG = timeG.right
  G = G.right
  errG = errG.right }

len = length(timeG)
limit2 = len   #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
for (i in 1:(len-1)){
  if(timeG[i+1]-timeG[i] >150) {
    limit2 = i  #make a new variable limit2 because limit already has a value
  }
}
timeG.left = timeG[1:limit2]
G.left = G[1:limit2]
errG.left = errG[1:limit2]
timeG.right = timeG[(limit2+1):len]
G.right = G[(limit2+1):len]
errG.right = errG[(limit2+1):len]
if (length(timeG.left)>length(timeG.right)){
  timeG = timeG.left
  G = G.left
  errG = errG.left }
if (length(timeG.right)>length(timeG.left)){
  timeG = timeG.right
  G = G.right
  errG = errG.right }


iqr = IQR(G)
sd = iqr*20/27
median = median(G)
upper = median + 3*sd
lower = median - 3*sd
timeG = timeG[G<upper & G>lower]
errG = errG[G<upper & G>lower]
G = G[G<upper & G>lower]


len = length(timeR)
limit = len    #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
for (i in 1:(len-1)){
  if(timeR[i+1]-timeR[i] >150) {
    limit = i
  }
}
timeR.left = timeR[1:limit]
R.left = R[1:limit]
errR.left = errR[1:limit]
timeR.right = timeR[(limit+1):len]
R.right = R[(limit+1):len]
errR.right = errR[(limit+1):len]
if (length(timeR.left)>length(timeR.right)){
  timeR = timeR.left
  R = R.left
  errR = errR.left }
if (length(timeR.right)>length(timeR.left)){
  timeR = timeR.right
  R = R.right
  errR = errR.right }

len = length(timeR)
limit2 = len   #in case the data is equally distributed through time (and therefore i is nan), this will help reserve the initial data set.
for (i in 1:(len-1)){
  if(timeR[i+1]-timeR[i] >150) {
    limit2 = i  #make a new variable limit2 because limit already has a value
  }
}
timeR.left = timeR[1:limit2]
R.left = R[1:limit2]
errR.left = errR[1:limit2]
timeR.right = timeR[(limit2+1):len]
R.right = R[(limit2+1):len]
errR.right = errR[(limit2+1):len]
if (length(timeR.left)>length(timeR.right)){
  timeR = timeR.left
  R = R.left
  errR = errR.left }
if (length(timeR.right)>length(timeR.left)){
  timeR = timeR.right
  R = R.right
  errR = errR.right }


iqr = IQR(R)
sd = iqr*20/27
median = median(R)
upper = median + 3*sd
lower = median - 3*sd
timeR = timeR[R<upper & R>lower]
errR = errR[R<upper & R>lower]
R = R[R<upper & R>lower]


period.G = 1/frequency
period.sec.G = period.G*24*60*60
phase1.G = (timeG%%period.G)/period.G

dataG = data.frame("mag"=G,"phase"=phase1.G,"err"=errG)
dataG = dataG[order(phase1.G),]
G = dataG$mag
phase1.G = dataG$phase
errG = dataG$err

phase2.G = phase1.G + 1
phase.G = c(phase1.G,phase2.G)
G.1 = G
#G.1 = G.1 - mean(G.1)
G = c(G,G)
#G = G - mean(G)


period.R = 1/frequency
period.sec.R = period.R*24*60*60
phase1.R = (timeR%%period.R)/period.R

dataR = data.frame("mag"=R,"phase"=phase1.R,"err"=errR)
dataR = dataR[order(phase1.R),]
R = dataR$mag
phase1.R = dataR$phase
errR = dataR$err

phase2.R = phase1.R + 1
phase.R = c(phase1.R,phase2.R)
R.1 = R
#R.1 = R.1 - mean(R.1)
R = c(R,R)
#R = R - mean(R)

ampl.G = max(G + errG) - min(G - errG)
ampl.R = max(R + errR) - min(R - errR)
ampl = max(ampl.G,ampl.R)
lower.G = mean(G) - alpha*ampl
upper.G = mean(G) + alpha*ampl
lower.R = mean(R) - alpha*ampl
upper.R = mean(R) + alpha*ampl

err = c(errG,errG,errR,errR)
errG = c(errG,errG)
errR = c(errR,errR)

df.G = data.frame("mag"=G,"phase"=phase.G,"band"="g")
df.R = data.frame("mag"=R,"phase"=phase.R,"band"="r")
df = rbind(df.G,df.R)

setwd(plot.dir)

filename.pc = paste("Phase Curve 2-band ",sourceid,".png",sep="")
png(file=filename.pc,width=700,height=520)
#par(mar=c(6,6,5,3)+.1)

title = paste("Gaia EDR3 ",sourceid,sep="")

if(period.sec.G<60){xlab = sprintf("phase (Period = %.3f secs)",period.sec.G)}
if(period.sec.G>=60 & period.sec.G<3600){xlab = sprintf("phase (Period = %.3f mins)",period.sec.G/60)}
if(period.sec.G>=3600 & period.sec.G<86400){xlab = sprintf("phase (Period = %.3f hours)",period.sec.G/3600)}
if(period.sec.G>=86400 & period.sec.G<31536000){xlab = sprintf("phase (Period = %.3f days)",period.sec.G/86400)}
if(period.sec.G>=31536000){xlab = sprintf("phase (Period = %.3f years)",period.sec.G/31536000)}

plotG = ggplot(df.G,aes(x=phase,y=mag)) + geom_point(size=2.5) + 
  #geom_line(linetype="dashed") + 
  theme_bw() +
  geom_errorbar(aes(ymin=G-errG,ymax=G+errG)) + ylim(upper.G,lower.G) + 
  ylab("ZTF g mag") +xlab(xlab) + 
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=27)) +
  theme(strip.text = element_text(size=20, face="bold")) + ggtitle(title) +
  theme(plot.title = element_text(size=27, face="bold",hjust=0.5)) +
  theme(plot.margin = margin(2, 55, 2, 45, "pt")) + 
  theme(legend.position = "none") + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x = element_blank()) 

plotR = ggplot(df.R,aes(x=phase,y=mag)) + geom_point(size=2.5,color="red") + 
  #geom_line(linetype="dashed",color="red") + 
  theme_bw() + ylim(upper.R,lower.R) +
  geom_errorbar(aes(ymin=R-errR,ymax=R+errR),color="red") + 
  ylab("ZTF r mag") +xlab(xlab) + 
  theme(axis.text=element_text(size=20,face="bold"),axis.title=element_text(size=27)) +
  theme(strip.text = element_text(size=20, face="bold")) +
  #theme(plot.title = element_text(size=27, face="bold")) +
  theme(plot.margin = margin(2, 55, 2, 45, "pt")) +
  theme(legend.position = "none") + theme(title = element_blank())

plot_grid(plotG, plotR, nrow = 2)


dev.off()


