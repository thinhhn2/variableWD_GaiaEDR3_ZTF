file = "C:/Work/Variable candidates/Variable_timeseries_g.csv"
file.query = "C:/Work/Variable candidates/Variable ZTF query g.csv"
fileR = "C:/Work/Variable candidates/Variable_timeseries_r.csv"
file.queryR = "C:/Work/Variable candidates/Variable ZTF query r.csv"

ztf.data = read.csv(file=file,header=TRUE,stringsAsFactor=TRUE)
ztf.query = read.csv(file=file.query,header=TRUE,stringsAsFactor=TRUE) 
ztf.dataR = read.csv(file=fileR,header=TRUE,stringsAsFactor=TRUE)
ztf.queryR= read.csv(file=file.queryR,header=TRUE,stringsAsFactor=TRUE) 

periodic.n = c(24,25,38,60,68,78,91,104,106,107,113,123,136,145,146,150,151,157,194,
               196,211,212,214,222,223,224,227,254,257,272,273,277,278,287,290,
               298,357,360,372,374,398,404,469,470,480,492,499,523,530,541,543,548,
               554,557,564,566,579,584,587,592,601,605,610,617,625,626,645,666,668,
               680,717,723,749,760,796,808,810,856,867,873,882,914,933,951,961,968,
               980,1008,1029,1044,1050,1051,1074,1080,1082,1098,1108,1110,1150,1154,1160,
               1168,1171,1174,1179,1209,1225,1228,1234,1287,1299,1313,1317,1318,
               1330,1334,1338,1343,1348,1349,1352,1356,1369,1371,1375,1381,1385,
               1395,1402,1406,1416,1429,1432,1438,32,65,221,228,250,316,325,345,452,477,
               497,526,551,581,773,811,854,875,990,1019,1026,1040,1046,1064,1176,
               1284,1339,1382,1394,1411,1414,1421,1439)

#levels = levels(factor(ztf.data$oid))
amplG.list = c()
amplR.list = c()

for (n in periodic.n){
  levels = ztf.query$oid[ztf.query$cntr_01 == n]
  levelsR = ztf.queryR$oid[ztf.queryR$cntr_01 == n]
  
  G = ztf.data$mag[ztf.data$oid == levels]
  catflag = ztf.data$catflag[ztf.data$oid == levels]
  sharp = ztf.data$sharp[ztf.data$oid == levels]
  limitmag = ztf.data$limitmag[ztf.data$oid == levels]
  
  G = G[catflag == 0 & abs(sharp) < 0.25 & G < limitmag - 1.0]
  
  R = ztf.dataR$mag[ztf.dataR$oid == levelsR]
  catflagR = ztf.dataR$catflag[ztf.dataR$oid == levelsR]
  sharpR = ztf.dataR$sharp[ztf.dataR$oid == levelsR]
  limitmagR = ztf.dataR$limitmag[ztf.dataR$oid == levelsR]
  
  R = R[catflagR == 0 & abs(sharpR) < 0.25 & R < limitmagR - 1.0]
  
  amplG = quantile(G,0.95) - quantile(G,0.05)
  amplG.list = c(amplG.list,amplG)
  
  amplR = quantile(R,0.95) - quantile(R,0.05)
  amplR.list = c(amplR.list,amplR)
}
