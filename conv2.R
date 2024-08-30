setwd("~/cwork/stl_tab/")
rm(list=ls())
scl = function(img,na.rm=TRUE) (img-min(img,na.rm=TRUE))/(max(img,na.rm=TRUE)-min(img,na.rm=TRUE))

case  = "sample_2"
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

X = cbind( 1, as.vector(row(tab2_s)*dx), as.vector(col(tab2_s)*dx), as.vector(tab2_s) )
dim(X)
#l = apply(X,2,function(x) sqrt(sum(x^2)))
#l = l/l[4]
#X = t(t(X)/l)

constr = matrix(c(0,0,0,1),1,4)
constr_rhs = 1
M = cbind(rbind(t(X) %*% X, constr),rbind(t(constr),0))

a1 = fft(tab1)
RHS = matrix(0,ncol(M),nrow(tab1)*ncol(tab1))
for (i in 1:4) {
  H = matrix(0,nrow(tab1),ncol(tab1))
  H[seq_len(nrow(tab2_s)),seq_len(ncol(tab2_s))] = X[,i]
  a2 = fft(H[c(1,nrow(H):2),c(1,ncol(H):2)])
  HC = Re(fft(a1*a2,inverse = TRUE))/length(a1)
  RHS[i,] = as.vector(HC)
}
RHS[seq_len(nrow(constr))+ncol(X),] = constr_rhs

BETA = solve(M, RHS)

a1 = fft(tab1^2)
H[seq_len(nrow(tab2_s)),seq_len(ncol(tab2_s))] = 1
a2 = fft(H[c(1,nrow(H):2),c(1,ncol(H):2)])
HC = Re(fft(a1*a2,inverse = TRUE))/length(a1)
Y2 = HC

rm(a1,a2,H,HC); gc()
R = Y2 - 2*colSums(RHS[1:4,] * (BETA[1:4,])) + colSums((t(X) %*% X %*% BETA[1:4,])*(BETA[1:4,]))

ret = R
ret[nrow(tab1) - seq_len(nrow(tab2_s))+1,] = NA
ret[,ncol(tab1) - seq_len(ncol(tab2_s))+1] = NA
i = which.min(as.vector(ret))
xshift = row(ret)[i] - 1
yshift = col(ret)[i] - 1
rm(ret); gc()
beta = BETA[,i]

img = -R
img = scl(img)
img[xshift+1+(-50:50),yshift+1] = 0
img[xshift+1,yshift+1+(-50:50)] = 0
png::writePNG(scl(img),paste0(case, "_shift.png"))

img = scl(tab1)
k=8
img[xshift+seq_len(nrow(tab2_s)),yshift+c(-k:k,-k:k+ncol(tab2_s))] = 0
img[xshift+c(-k:k,-k:k+nrow(tab2_s)),yshift+seq_len(ncol(tab2_s))] = 0
png::writePNG(scl(img),paste0(case1, ".png"))

tab1_sel = tab1[seq_len(nrow(tab2_s))+xshift,seq_len(ncol(tab2_s))+yshift]
img = tab1_sel
png::writePNG(scl(img),paste0(case1, "_cut.png"))
img = tab2
png::writePNG(scl(img),paste0(case2, ".png"))

img = tab1_sel - tab2_s
png::writePNG(scl(img),paste0(case1, "_diff.png"))

img[] = X %*% beta[1:4]
png::writePNG(scl(img),paste0(case2, "_even.png"))

i = ncol(tab1_sel)*0.5
matplot(cbind(img[,i], tab1_sel[,i], tab2_s[,i]+beta[1]),type="l",lty=1)

i = seq_len(nrow(tab2))
j = ncol(tab2)/2
x = i*dx/4
y = j*dx/4
z2 = tab2[i,j]
az = beta[1]+beta[2]*x+beta[3]*y + z2
plot(x,az,type="l")
i = unique(i/4)
j = unique(j/4)
x = i*dx
z1 = tab1[i+xshift,j+yshift]
lines(x,z1,col=3)
