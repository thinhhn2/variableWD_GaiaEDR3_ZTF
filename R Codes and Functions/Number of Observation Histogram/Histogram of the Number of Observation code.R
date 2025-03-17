library(ggplot2)

gaia.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Data 2000 and 300 Targets/Gaia EDR3 Data/2198 candidates.csv"
ztf.g.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Data 2000 and 300 Targets/ZTF DR5/ZTF DR5 g number of observations nobsrel.csv"
ztf.r.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Data 2000 and 300 Targets/ZTF DR5/ZTF DR5 r number of observations nobsrel.csv"
plot.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Data 2000 and 300 Targets/Number of observations plot"

ztf.g = read.csv(ztf.g.file,header=TRUE)
ztf.r = read.csv(ztf.r.file,header=TRUE)
gaia = read.csv(gaia.file,header=TRUE,colClasses=c("source_id"="character"))

#Import Gaia number of observations
conver_factor_g = 8.86
G = gaia$phot_g_n_obs/conver_factor_g
BP = gaia$phot_bp_n_obs
RP = gaia$phot_rp_n_obs

#Import ZTF number of observations
g = ztf.g$nobsrel
r = ztf.r$nobsrel

#Convert histogram to bar plot data
start.x = 0
end.x = (max(c(g,r))%/%100 + 1)*100
by.x = 10
xrange = seq(from=start.x,to=end.x,by=by.x)
G.bp = hist(G,xrange)$counts
BP.bp = hist(BP,xrange)$counts
RP.bp = hist(RP,xrange)$counts
g.bp = hist(g,xrange)$counts
r.bp = hist(r,xrange)$counts


#Merge the data after numobs = 500 into one bin
lim = 500
lim.index = 500/by.x + 1
  
for (i in (lim.index+1):length(G.bp)){ G.bp[lim.index] = G.bp[lim.index] + G.bp[i]}
G.bp = G.bp[c(1:lim.index)]
for (i in (lim.index+1):length(BP.bp)){ BP.bp[lim.index] = BP.bp[lim.index] + BP.bp[i]}
BP.bp = BP.bp[c(1:lim.index)]
for (i in (lim.index+1):length(RP.bp)){ RP.bp[lim.index] = RP.bp[lim.index] + RP.bp[i]}
RP.bp = RP.bp[c(1:lim.index)]
for (i in (lim.index+1):length(g.bp)){ g.bp[lim.index] = g.bp[lim.index] + g.bp[i]}
g.bp = g.bp[c(1:lim.index)]
for (i in (lim.index+1):length(r.bp)){ r.bp[lim.index] = r.bp[lim.index] + r.bp[i]}
r.bp = r.bp[c(1:lim.index)]

#Make histogram 
xaxis = xrange+5
xaxis = xaxis[c(1:lim.index)]
df.G = data.frame("numobs"=xaxis,"counts"=G.bp,"band"="Gaia EDR3 G")
df.BP = data.frame("numobs"=xaxis,"counts"=BP.bp,"band"="Gaia EDR3 BP")
df.RP = data.frame("numobs"=xaxis,"counts"=RP.bp,"band"="Gaia EDR3 RP")
df.g = data.frame("numobs"=xaxis,"counts"=g.bp,"band"="ZTF DR5 g")
df.r = data.frame("numobs"=xaxis,"counts"=r.bp,"band"="ZTF DR5 r")

df = rbind(df.RP,df.r,df.BP,df.g,df.G)

setwd(plot.dir)
filename = sprintf("Histogram of the Number of Observation from Gaia and ZTF.png")
png(file=filename,width=750,height=700)

plot = ggplot(data=df, aes(x=numobs, y=counts, fill=band, color=band)) + geom_bar(stat="identity") + facet_wrap(vars(band),ncol=2) + scale_fill_manual(values=c("light coral","red","light blue","blue","forest green")) + scale_color_manual(values=c("black", "black", "black","black","black")) + theme_bw() +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) + theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=22,face="bold")) + 
  theme(legend.text = element_text(size=18,face="bold"),legend.title = element_blank()) + xlab("Number of Observations") + 
  theme(legend.position="top") + theme(strip.text = element_text(size=20, face="bold"))

#ggsave(file="Histogram of the Number of Observation from Gaia and ZTF.svg",plot=plot,width=11,height=8)
dev.off()
