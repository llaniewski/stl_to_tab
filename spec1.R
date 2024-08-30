library(rfracture)

tab1_b = matrix(0,nrow(tab1)*2,ncol(tab1)*2)
tab1_b[seq_len(nrow(tab1)),seq_len(ncol(tab1))] = tab1
a = fft(tab1_b)
rm(tab1_b); gc()
tab = expand.grid(fx=seq_circ(nrow(a)),fy=seq_circ(ncol(a)))
tab$val = as.vector(Re(a))
rm(a); gc()

tab$f = sqrt(tab$fx^2+tab$fy^2)
range(tab$f)
tab$cf = cut(tab$f,breaks=seq(0,8000,100))

ret = by(seq_len(nrow(tab)), tab$cf, sample, 100)
ret = do.call(c,ret)
tab_samp = tab[ret,]

dim(tab_samp)
plot(tab_samp$f, tab_samp$val,log="x")

library(gamlss)
mod = gamlss(val~0, sigma.formula = ~log(f),data=tab_samp)
plot(mod)
