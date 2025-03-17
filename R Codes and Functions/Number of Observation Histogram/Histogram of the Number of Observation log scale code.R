library(ggplot2)

gaia.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Uniform Distribution 2423 targets/Gaia data 2423 targets.csv"
ztf.g.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Number of observations plot/ZTF 2423 query g band.csv"
ztf.r.file = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Number of observations plot/ZTF 2423 query r band.csv"
plot.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Comparison between Gaia,ZTF, and CRTS/Number of observations plot"

ztf.g = read.csv(ztf.g.file,header=TRUE,stringsAsFactor=FALSE)
ztf.r = read.csv(ztf.r.file,header=TRUE,stringsAsFactor=FALSE)
gaia = read.csv(gaia.file,header=TRUE,colClasses=c("source_id"="character"))

#Import Gaia number of observations
conver_factor_g = 8.86
G = gaia$phot_g_n_obs/conver_factor_g
BP = gaia$phot_bp_n_obs
RP = gaia$phot_rp_n_obs

#Import ZTF number of observations
g = as.numeric(ztf.g$nobsrel[ztf.g$nobsrel != "null"])
r = as.numeric(ztf.r$nobsrel[ztf.r$nobsrel != "null"])


df.G = data.frame("counts"=G,"band"="Gaia DR3 G")
df.BP = data.frame("counts"=BP,"band"="Gaia DR3 BP")
df.RP = data.frame("counts"=RP,"band"="Gaia DR3 RP")
df.g = data.frame("counts"=g,"band"="ZTF DR8 g")
df.r = data.frame("counts"=r,"band"="ZTF DR8 r")

df = rbind(df.RP,df.r,df.BP,df.g,df.G)

setwd(plot.dir)
filename = sprintf("Histogram of the Number of Observation from Gaia DR3 and ZTF DR8 log scale.png")
png(file=filename,width=750,height=700)

ggplot(data=df, aes(x=counts, fill=band, color=band)) + geom_histogram(position="identity") + facet_wrap(vars(band),ncol=2) + scale_x_log10() + scale_fill_manual(values=c("red","light coral","blue","light blue","forest green")) + scale_color_manual(values=c("black", "black", "black","black","black")) + theme_bw() +
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(axis.text=element_text(size=16,face="bold"),axis.title=element_text(size=22,face="bold")) + 
  #theme(legend.text = element_text(size=18,face="bold"),legend.title = element_blank()) +
  xlab("Number of Observations") + 
  theme(legend.position="none") + theme(strip.text = element_text(size=20, face="bold"))

#ggsave(file="Histogram of the Number of Observation from Gaia and ZTF.svg",plot=plot,width=11,height=8)
dev.off()
