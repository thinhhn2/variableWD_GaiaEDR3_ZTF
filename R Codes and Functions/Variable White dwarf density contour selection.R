data = read.csv(file=file.choose(),colClasses=c("source_id"="character"))
library("MASS")
library("sp")

bprp = data$phot_bp_mean_mag - data$phot_rp_mean_mag
absG = data$phot_g_mean_mag + 5*log10(data$parallax/1000) + 5
dens = kde2d(bprp,absG,n=500)

plot(bprp,absG,ylim=c(max(absG),min(absG)),pch=16)
prob.list = c(0.1,0.2,0.3,0.4,0.5,0.6,0.65,0.68,0.7,0.8,0.9)
#prob.list = c(0.48,0.5)
color.list = 2:12
for (i in 1:length(prob.list)){
prob = prob.list[i]
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})
contour(dens, levels=levels, labels=prob, add=T, col=color.list[i], lwd =2)
}

#Probability contour density for ZZ Cet
prob = 0.1
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels2 <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})

ls <- contourLines(dens, level=levels2)
zzcet <- point.in.polygon(bprp, absG, ls[[1]]$x, ls[[1]]$y)

zzcet_candidates = data[zzcet==1,]


data = read.csv(file=file.choose(),colClasses=c("source_id"="character"))
library("MASS")
library("sp")

bprp = data$phot_bp_mean_mag - data$phot_rp_mean_mag
absG = data$phot_g_mean_mag + 5*log10(data$parallax/1000) + 5
dens = kde2d(bprp,absG,n=500)

plot(bprp,absG,ylim=c(max(absG),min(absG)),pch=16)
prob.list = c(0.1,0.15,0.2,0.3,0.4,0.45,0.5,0.6,0.7,0.8,0.9)
#prob.list = c(0.48,0.5)
color.list = 2:12
for (i in 1:length(prob.list)){
prob = prob.list[i]
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})
contour(dens, levels=levels, labels=prob, add=T, col=color.list[i], lwd =2)
}


#Probability contour density for V777 Herculis
prob = 0.45
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels1 <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})

ls <- contourLines(dens, level=levels1)
v777her <- point.in.polygon(bprp, absG, ls[[2]]$x, ls[[2]]$y)
v777her_candidates = data[v777her==1,]

#Probability contour density for GW Virginis
prob = 0.45
dx <- diff(dens$x[1:2])
dy <- diff(dens$y[1:2])
sz <- sort(dens$z)
c1 <- cumsum(sz) * dx * dy 
levels1 <- sapply(prob, function(x) { 
    approx(c1, sz, xout = 1 - x)$y
})

ls <- contourLines(dens, level=levels1)
gwvir <- point.in.polygon(bprp, absG, ls[[1]]$x, ls[[1]]$y)

gwvir_candidates = data[gwvir==1,]

data = read.csv(file=file.choose(),colClasses=c("source_id"="character"))
bprp = data$phot_bp_mean_mag - data$phot_rp_mean_mag
absG = data$phot_g_mean_mag + 5*log10(data$parallax/1000) + 5

plot(bprp,absG,ylim=c(max(absG),min(absG)),pch=16)
points(zzcet_candidates$phot_bp_mean_mag - zzcet_candidates$phot_rp_mean_mag,zzcet_candidates$phot_g_mean_mag + 5*log10(zzcet_candidates$parallax/1000) + 5,col="red",pch=16)
points(v777her_candidates$phot_bp_mean_mag - v777her_candidates$phot_rp_mean_mag,v777her_candidates$phot_g_mean_mag + 5*log10(v777her_candidates$parallax/1000) + 5,col="blue",pch=16)
points(gwvir_candidates$phot_bp_mean_mag - gwvir_candidates$phot_rp_mean_mag,gwvir_candidates$phot_g_mean_mag + 5*log10(gwvir_candidates$parallax/1000) + 5,col="green",pch=16)


