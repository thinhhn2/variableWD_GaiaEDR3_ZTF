ztf.pvalue = function(file){
#oid.list = c()
ra.list = c()
dec.list = c()
meanmag.list = c()
pvalue.list = c()

ztf.data = read.csv(file=file,header=TRUE,stringsAsFactor=TRUE)
levels = levels(factor(ztf.data$oid))

for (i in 1:length(levels)) {
mag = ztf.data$mag[ztf.data$oid == levels[i]]
time = ztf.data$mjd[ztf.data$oid == levels[i]]
err = ztf.data$magerr[ztf.data$oid == levels[i]]
catflag = ztf.data$catflag[ztf.data$oid == levels[i]]
sharp = ztf.data$sharp[ztf.data$oid == levels[i]]
limitmag = ztf.data$limitmag[ztf.data$oid == levels[i]]
ra = ztf.data$ra[ztf.data$oid == levels[i]]
dec = ztf.data$dec[ztf.data$oid == levels[i]]

ra = ra[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
dec = dec[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
time = time[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
err = err[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]
mag = mag[catflag == 0 & abs(sharp) < 0.25 & mag < limitmag - 1.0]

#err = 0*err + 0.01
#mag = rnorm(length(err),mean=0,sd=0.01)

mean.ra = mean(ra)
mean.dec = mean(dec)
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
ra.list = c(ra.list,mean.ra)
dec.list = c(dec.list,mean.dec)
}

return(data.frame("ra"=ra.list,"dec"=dec.list,"meanmag"=meanmag.list,"pvalue"=pvalue.list))
}
