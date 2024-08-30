rm(list=ls())

source("lib.R")

froll = 5
p.iso = function(f) ifelse(f<froll, 1, (f/froll)^{-3})
span = c(4,4)
frac = fracture_matrix(dims = span*100, span = diag(span), period = diag(span), power.iso = p.iso)
tab = frac$f1
plot(frac)

x = seq(0,4,len=100)
plot(x,sin(((x/4-0.5)^3*pi)))
plot(x,pnorm((x/4-0.5)*10)-x/4)

x = (frac$points[,1]-2)+(frac$points[,2]-2)*pi
d = pnorm(x*10)-x/4

tab = frac$f1+d*100
plot(tab[1,])


dim(tab)
tab = tab[1:100+100,1:200+200]
span = diag(frac$span)/dim(frac$f1) * dim(tab)

spec = spec2d(tab, span = span, plot=FALSE)

i = binned.sample(spec$f)
plot(spec$f[i],spec$p[i],log="xy",pch=16,cex=0.4)
f = seq(0,300,len=100)
lines(f,p.iso(f),col=2)



