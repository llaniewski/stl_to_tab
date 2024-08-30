

for (case2 in c("sample_1_mag10","sample_1_mag40","sample_2_mag10","sample_2_mag40","sample_3_mag10")) {
load(paste0("~/cwork/stl_tab/",case2,".Rdata"))
dx = mean(attr(tab,"pixe"))
prof = tab[,ncol(tab)/2]
png(paste0(case2,"_prof1.png"),width=600, height=600, res = 150)
plot(1e3*seq_len(nrow(tab))*dx, 1e3*prof,asp=10,type="l",xlab="X [mm]",ylab="Y [mm]",main="Profile")
dev.off()
png(paste0(case2,"_prof1_spec.png"),width=600, height=600, res = 150)
prof = ts(prof,deltat=dx)
sp = spectrum(prof,col=8)
m = lm(log(sp$spec)~log(sp$freq))
x = sp$freq
y = x^(m$coefficients[2])*exp(m$coefficients[1])
lines(x,y,col=1,lty=2)
text(x[length(x)*0.7],y[length(x)*0.7],sprintf("~ f^(%.2f)",m$coefficients[2]),adj=c(0.2,-1))
dev.off()
}

sp = spectrum(prof,col=8)



case2

case2 = "sample_1_mag10"
prof = 1e3*tab[,floor(ncol(tab)*seq(0,1,len=20))]
matplot(1e3*seq_len(nrow(tab))*dx, t(t(prof)+0.3*(-1:1)),asp=10,type="l",xlab="X [mm]",ylab="Y [mm]",lty=1,col=1)
prof = ts(prof,deltat=dx)
sp = spectrum(prof,col=8,lty=1)

