setwd("~/cwork/stl_tab/")
rm(list=ls())
scl = function(img) (img-min(img))/(max(img)-min(img))

case  = "sample_1"
case1 = paste0(case, "_mag10")
case2 = paste0(case, "_mag40")

load(paste0("~/cwork/stl_tab/",case1,".Rdata"))
tab1 = tab
load(paste0("~/cwork/stl_tab/",case2,".Rdata"))
tab2 = tab
rm(tab)

dx = mean(attr(tab1,"pixel"))
subsamp = round(mean(attributes(tab1)$pixel/attributes(tab2)$pixel))

tab2_s = tab2[seq(1,nrow(tab2),subsamp),seq(1,ncol(tab2),subsamp)]

tab2_sb = matrix(0,nrow(tab1),ncol(tab1))
tab2_sb[seq_len(nrow(tab2_s)),seq_len(ncol(tab2_s))] = tab2_s
tab2_1 = matrix(0,nrow(tab1),ncol(tab1))
tab2_1[seq_len(nrow(tab2_s)),seq_len(ncol(tab2_s))] = 1

a1 = fft(tab1)
a2 = fft(tab2_sb[c(1,nrow(tab2_sb):2),c(1,ncol(tab2_sb):2)])
a3 = fft(tab2_1[c(1,nrow(tab2_1):2),c(1,ncol(tab2_1):2)])
ret_sb   = fft(a1*a2,inverse = TRUE)/length(a1)
ret_1 = fft(a1*a3,inverse = TRUE)/length(a1)

ret_sb = ret_sb[seq_len(nrow(tab1)-nrow(tab2_s)),seq_len(ncol(tab1)-ncol(tab2_s))]
ret_1 = ret_1[seq_len(nrow(tab1)-nrow(tab2_s)),seq_len(ncol(tab1)-ncol(tab2_s))]

ret = ret_sb - ret_1*sum(tab2_sb)/sum(tab2_1)

img = Re(ret)

i = which.max(as.vector(img))
xshift = row(ret)[i] - 1
yshift = col(ret)[i] - 1

img = scl(img)
img[xshift+1+(-50:50),yshift+1] = 0
img[xshift+1,yshift+1+(-50:50)] = 0
png::writePNG(scl(img),paste0(case, "_shift.png"))

img = scl(tab1)
k=8
img[xshift+seq_len(nrow(tab2_s)),yshift+c(-k:k,-k:k+ncol(tab2_s))] = 0
img[xshift+c(-k:k,-k:k+nrow(tab2_s)),yshift+seq_len(ncol(tab2_s))] = 0
png::writePNG(scl(img),paste0(case1, ".png"))

tab2_sbs = tab2_sb[
  (seq_len(nrow(tab2_sb))-xshift-1) %% nrow(tab2_sb) + 1,
  (seq_len(ncol(tab2_sb))-yshift-1) %% ncol(tab2_sb) + 1]
img = tab2_sbs
png::writePNG(scl(img),paste0(case2, "_big.png"))

tab1_sel = tab1[seq_len(nrow(tab2_s))+xshift,seq_len(ncol(tab2_s))+yshift]
img = tab1_sel
png::writePNG(scl(img),paste0(case1, "_cut.png"))
img = tab2
png::writePNG(scl(img),paste0(case2, ".png"))

img = tab1_sel - tab2_s
png::writePNG(scl(img),"test_s_diff.png")

tab = data.frame(
  x  = as.vector(dx * row(tab2_s)),
  y  = as.vector(dx * col(tab2_s)),
  v1 = as.vector(tab1_sel),
  v2 = as.vector(tab2_s)
)
m = lm(I(v1-v2) ~ x+y,data=tab)
summary(m)

tab2_adj = tab2_s + matrix(m$fitted.values, nrow(tab2_s), ncol(tab2_s))
png::writePNG(scl(img),"test_s2_adj.png")

tab2_adj = tab2 + m$coefficients[1] + m$coefficients[2] * (dx/subsamp) * row(tab2) + m$coefficients[3] * (dx/subsamp) * col(tab2)
img = tab2_adj 
png::writePNG(scl(img),"test_s2_adj.png")

i = 70

plot(seq_len(nrow(tab1))*dx, tab1[,i+yshift],asp=10,type="l",xlim=c(xshift*dx,xshift*dx+nrow(tab2)*dx/subsamp))
lines((seq_len(nrow(tab2)))*dx/subsamp+xshift*dx, tab2_adj[,i*subsamp],col=3)

i = 30
plot(seq_len(nrow(tab1))*attr(tab1,"pixel")[1], tab1[,i+yshift],asp=50,type="l")
lines((seq_len(nrow(tab2)))*attr(tab2,"pixel")[1]+xshift*attr(tab1,"pixel")[1], tab2[,4*i]+zshift,col=4)


plot(tab2[,i*4]+zshift)
lines(tab1_sel[,i])

plot(tab2_s[,30]+zshift)
lines(tab1_sel[,30])

mean(tab2_s+zshift - tab1_sel)


dim(tab1)*dx
dim(tab2)*dx/4
