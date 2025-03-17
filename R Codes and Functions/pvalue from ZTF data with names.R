#data.dir = "C:/Thinh Nguyen files/Villanova University/Summer Research 2020 Switzerland/Gaia EDR3 Data/Gaia EDR3 Variable White Dwarf Candidates Reselection/Reselection with ZTF DR5/V777 Her/Data"

ztf.pvalue = function(name,band,data.dir){

library(FITSio)
band.list = c()
meanmag.list = c()
pvalue.list = c()

for (n in name){

file = sprintf("%s%s.fits",n,band)

setwd(data.dir)

Data = readFITS(file)
mag = Data$col[[6]]
time = Data$col[[5]]
err = Data$col[[7]]
catflag = Data$col[[8]]
sharp = Data$col[[13]]
limitmag = Data$col[[18]]

time = time[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
err = err[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
mag = mag[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]

#err = 0*err + 0.01
#mag = rnorm(length(err),mean=0,sd=0.01)

mean = weighted.mean(mag,(1/err^2))
df = length(mag) - 1
cq.stat = 0

for (i in 1:length(mag)){
cq.stat = cq.stat + ((mag[i] - mean)^2/err[i]^2)
}

p.value = pchisq(cq.stat,df,lower.tail=FALSE)

if (is.nan(mean) == TRUE) {p.value = NA}

#oid.list = c(oid.list,levels[i])
meanmag.list = c(meanmag.list,mean)
pvalue.list = c(pvalue.list,p.value)
}

return(data.frame("name"=name,"meanmag"=meanmag.list,"pvalue"=pvalue.list))
}
