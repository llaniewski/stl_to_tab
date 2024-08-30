library(raster)
library(rfracture)

setwd("~/cwork/stl_tab/")
rm(list=ls())
source("lib.R")


case  = "sample_2"
case1 = paste0(case, "_mag10")
case2 = paste0(case, "_mag40")

load(paste0("~/cwork/stl_tab/",case1,".Rdata"))
tab1 = tab
load(paste0("~/cwork/stl_tab/",case2,".Rdata"))
#tab2 = tab[,1:floor(ncol(tab)/2)]
#attr(tab2,"pixel") = attr(tab,"pixel")
tab2 = tab
rm(tab)

scale = 1e3
tab1 = tab1*scale
span1 = dim(tab1) * attr(tab1,"pixel") * scale
tab2 = tab2*scale
span2 = dim(tab2) * attr(tab2,"pixel") * scale

spec = spec2d(tab1,span1,plot=FALSE)
i = binned.sample(spec$f)
spec1 = spec[i,]
spec = spec2d(tab2,span2,plot=FALSE)
i = binned.sample(spec$f)
spec2 = spec[i,]

#spec = spec2d(tab1,span1)
flim1 = 0.25*nrow(tab1)/(span1[1])
flim2 = 0.25*nrow(tab2)/(span2[1])
sp = rbind(cbind(spec1[spec1$f<flim1,],tab="tab1"),cbind(spec2[spec2$f<flim2,],tab="tab2"))
sp$tab = factor(sp$tab)
m = lm(log(p)~log(f)+tab,data=sp)

f = seq(0,max(spec$f),len=200)
plot(spec2$f,spec2$p,log="xy",type="n",xlab="Freq [1/mm]", ylab="2D Power Spectrum [mm^4]")
points(spec1$f,spec1$p,pch=16,cex=0.2,col=2)
points(spec2$f,spec2$p,pch=16,cex=0.2,col=3)
abline(m$coefficients[1]/log(10),m$coefficients[2])
abline((m$coefficients[1] + m$coefficients[3])/log(10),m$coefficients[2])

nms = c("x10 magnification", "x40 magnification")
png(paste0(case,"_spec2d.png"), width = 18, height = 6,units = "cm", res=300,pointsize = 8)
par(mfrow=c(1,3),mar=c(5.1,4.1,4.1,0.1))
tab_r = raster(tab1, xmn=0, xmx=span1[2], ymn=0, ymx=span1[1])
plot(tab_r,xlab="X [mm]", ylab="Y [mm]", legend = FALSE, main=nms[1])
tab_r = raster(tab2, xmn=0, xmx=span2[2], ymn=0, ymx=span2[1])
plot(tab_r,xlab="X [mm]", ylab="Y [mm]", legend = FALSE, main=nms[2])
plot(spec2$f,spec2$p,log="xy",type="n",xlab="Freq [1/mm]", ylab="2D Power Spectrum [mm^4]")
points(spec1$f,spec1$p,pch=16,cex=0.2,col=2)
points(spec2$f,spec2$p,pch=16,cex=0.2,col=3)
#abline(m$coefficients[1]/log(10),m$coefficients[2],lty=2)
#abline((m$coefficients[1] + m$coefficients[3])/log(10),m$coefficients[2],lty=2)
abline(coef[1]/log(10),coef[2],lty=2)
#legend("bottomleft", c(nms,sprintf("f^%.2f",m$coefficients[2])), pch=c(16,16,NA),col=c(2,3,1),lty=c(NA,NA,2))
legend("bottomleft", c(nms,sprintf("f^%.2f",coef[2])), pch=c(16,16,NA),col=c(2,3,1),lty=c(NA,NA,2))
dev.off()

####

coef = m$coefficients
power.iso = function(f)  exp(coef[1])*f^{coef[2]}
n = 250
span = mean(span1)
obj = fracture_matrix(dims = c(n,n), span = diag(2)*span, period = diag(2)*span, power.iso = power.iso, seed = 123)
tab1 = obj$f1
span1 = c(span,span)
obj = fracture_matrix(dims = c(n,n)*4, span = diag(2)*span, period = diag(2)*span, power.iso = power.iso, seed = 123)
tab2 = obj$f1[seq_len(n),seq_len(n)]
span2 = c(span,span)/4
case = paste0(case,"synth")


tab1 = obj$f1

tab2 = obj$f1[seq_len(n/4),seq_len(n/8)]

obj_r = raster(obj$f1[seq_len(n/4),seq_len(n/4)], xmn=0, xmx=span/4, ymn=0, ymx=span/4)
plot(obj_r,xlab="X [mm]", ylab="Y [mm]", legend = FALSE)


obj = fracture_matrix(dims = c(n,n), span = diag(2)*span/4, period = diag(2)*span, power.iso = power.iso, seed = 123)
obj_r = raster(obj$f1, xmn=0, xmx=obj$span[2,2], ymn=0, ymx=obj$span[1,1])
plot(obj_r,xlab="X [mm]", ylab="Y [mm]", legend = FALSE)





plot(seq_len(ncol(tab1))/ncol(tab1)*span1[2],tab1[nrow(tab1)/2,] - mean(tab1[nrow(tab1)/2,]),type="l")
lines(seq_len(ncol(obj$f1))/ncol(obj$f1)*span1[2], obj$f1[nrow(obj$f1)/2,])

########

val   = ts(tab1[,floor(ncol(tab1)/2)], deltat = attr(tab1,"pixel")[1])
spec1 = spectrum(val,plot=FALSE)
val   = ts(tab2[,floor(ncol(tab2)/2)], deltat = attr(tab2,"pixel")[1])
spec2 = spectrum(val,plot=FALSE)
plot(spec2$freq,spec2$spec,log="xy",col=2,pch=16,cex=0.3)
points(spec1$freq,spec1$spec,log="xy",col=3,pch=16,cex=0.3)

plot(tab1[floor(nrow(tab1)/2),])
plot(tab2[floor(nrow(tab1)/2),])
plot(val)

#####


flims = range(spec$f[-1])
#flims[2] = 0.5*nrow(tab2)/span2[1]
fc = cut(log(spec$f),breaks = seq(log(flims[1]),log(flims[2]),len=30))
k = 500
i = by(seq_len(nrow(spec)), fc, function(x) if (length(x) < k) x else sample(x, size=k))
i = sort(do.call(c,i))

i = binned.sample(spec$f)
plot(spec$f[i], spec$p[i],log="xy",pch=16,cex=0.5)


#####

val = ts(tab1[1,],deltat = attr(tab1,"pixel")[2]*1e3)
spectrum(val)
val = ts(tab2[1,],deltat = attr(tab2,"pixel")[2]*1e3)
spectrum(val)



val   = ts(tab1[1,], deltat = attr(tab1,"pixel")[2])
spec1 = spectrum(val,plot=FALSE)
val   = ts(tab2[1,], deltat = attr(tab2,"pixel")[2])
spec2 = spectrum(val,plot=FALSE)


plot(spec2)
lines(spec1$freq,spec1$spec,col=2)

ret1 = lapply(floor(seq(1,nrow(tab1),len=500)), function(i) {
  val   = ts(tab1[i,], deltat = attr(tab1,"pixel")[2])
  spectrum(val,plot=FALSE)
})
spec1 = rowMeans(sapply(ret1,function(x) x$spec))
freq1 = rowMeans(sapply(ret1,function(x) x$freq))

ret2 = lapply(floor(seq(1,nrow(tab2),len=500)), function(i) {
  val   = ts(tab2[i,], deltat = attr(tab2,"pixel")[2])
  spectrum(val,plot=FALSE)
})
spec2 = rowMeans(sapply(ret2,function(x) x$spec))
freq2 = rowMeans(sapply(ret2,function(x) x$freq))

plot(freq2,spec2,log="xy",type="l",ylim=c(1e-23,1e-11))
lines(freq1,spec1,col=2)

y = tab1[1,]
dx = attr(tab1,"pixel")[2]
val   = ts(y, deltat = dx)
spec1 = spectrum(val,plot=FALSE)
plot(spec1)
y2 = spec.taper(y,p=0.1)
p = fft(y2)
p = (Re(p)^2+Im(p)^2)*dx/length(y)
f = rfracture::seq_circ(length(y)) / (length(y)*dx)
plot(spec1)
lines(abs(f),p,col=3)


plot(y)
plot(y2)

spec.taper

library(rfracture)
n = length(y)
taper.w = function(n,p=0.1) {
  w = function(w) {ifelse(abs(w)<1,0.5*(1-cos(w*pi)),1)}
  w(seq_circ(n)/(n*p))
}


wx = taper.w(nrow(tab1))
wy = taper.w(ncol(tab1))

tab1_taper = tab1 * outer(wx,wy)

tab1_p = fft(tab1_taper)
tab1_p2 = Re(tab1_p*Conj(tab1_p))*(dx/length(y))
tab1_f = sqrt(outer(seq_circ(nrow(tab1))^2,seq_circ(ncol(tab1))^2,"+"))/ (length(y)*dx)
tab1_sp = data.frame( p = as.vector(tab1_p2), f = as.vector(tab1_f) )

i = sample.int(nrow(tab1_sp),10000)
plot(tab1_sp$f[i],tab1_sp$p[i],log="y")



