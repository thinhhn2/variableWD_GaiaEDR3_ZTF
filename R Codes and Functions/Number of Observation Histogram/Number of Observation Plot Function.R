library(ggplot2)
xrange = seq(from=0,to=900,by=10)
crts.bp = hist(crts[,1],xrange)$counts
gaia.g.bp = hist(gaia.g[,1],xrange)$counts
gaia.rp.bp = hist(gaia.rp[,1],xrange)$counts
gaia.bp.bp = hist(gaia.bp[,1],xrange)$counts
ztf.g.bp = hist(ztf.g[,1],xrange)$counts
ztf.r.bp = hist(ztf.r[,1],xrange)$counts

for (i in 51:90){ crts.bp[50] = crts.bp[50] + crts.bp[i]}
crts.bp = crts.bp[c(1:50)]
gaia.g.bp = gaia.g.bp[c(1:50)]
gaia.rp.bp = gaia.rp.bp[c(1:50)]
gaia.bp.bp = gaia.bp.bp[c(1:50)]
for (i in 51:90){ ztf.g.bp[50] = ztf.g.bp[50] + ztf.g.bp[i]}
ztf.g.bp = ztf.g.bp[c(1:50)]
for (i in 51:90){ ztf.r.bp[50] = ztf.r.bp[50] + ztf.r.bp[i]}
ztf.r.bp = ztf.r.bp[c(1:50)]

xaxis = xrange+5
xaxis = xaxis[c(1:50)]
df.crts = data.frame("numobs"=xaxis,"counts"=crts.bp,"band"="CRTS")
df.gaia.g = data.frame("numobs"=xaxis,"counts"=gaia.g.bp,"band"="Gaia G")
df.gaia.bp = data.frame("numobs"=xaxis,"counts"=gaia.bp.bp,"band"="Gaia BP")
df.gaia.rp = data.frame("numobs"=xaxis,"counts"=gaia.rp.bp,"band"="Gaia RP")
df.ztf.g = data.frame("numobs"=xaxis,"counts"=ztf.g.bp,"band"="ZTF g")
df.ztf.r = data.frame("numobs"=xaxis,"counts"=ztf.r.bp,"band"="ZTF r")

df = rbind(df.gaia.g,df.ztf.g,df.gaia.rp,df.ztf.r,df.gaia.bp,df.crts)


filename = sprintf("Histogram of the Number of Observation from 3 Surveys ver 3.png")
png(file=filename,width=700,height=700)


plot = ggplot(data=df, aes(x=numobs, y=counts, fill=band, color=band)) + geom_bar(stat="identity") + facet_wrap(vars(band),ncol=2) + scale_fill_manual(values=c("forest green","lightgreen","red","light coral","blue", "yellow")) + scale_color_manual(values=c("black", "black", "black","black","black","black")) + theme_bw() +
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +  ggtitle("Histogram of the \nNumber of Observation from the 3 Surveys") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(face = "bold",size=24)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=20,face="bold")) + 
	theme(legend.text = element_text(size=16,face="bold"),legend.title = element_text(size=16,face="bold")) + xlab("Number of Observation") + 
	theme(legend.position="top") + theme(strip.text = element_text(size=16, face="bold"))


par(mar=c(6,6,5,3)+.1)

crts.plot = ggplot(data=df.crts, aes(x=numobs, y=counts, fill=band, color=band)) + geom_bar(stat="identity") + scale_fill_manual(values="lightgreen") + scale_color_manual(values="black") 
gaia.plot = ggplot(data=df.gaia, aes(x=numobs, y=counts, fill=band, color=band)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("green", "blue", "red")) + scale_color_manual(values=c("black", "black", "black"))
ztf.plot = ggplot(data=df.ztf, aes(x=numobs, y=counts, fill=band, color=band)) + geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c("darkgreen", "darkred")) + scale_color_manual(values=c("black", "black"))

main = "Histogram of the \nObservation Times from the 3 Surveys"
barplot(crts.bp,col="light coral",names.arg=xrange[-1],space=0,las=1,main=main,cex.main=2.2,cex.names=1.8,cex.axis=1.8)
title(ylab="Counts",cex.lab=2,line=4)
legend("topleft",legend="CRTS DR2",pch=16,lty=2,col="light coral",cex=1.6)
dev.off()