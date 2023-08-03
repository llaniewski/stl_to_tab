
library(Rcpp)

Rcpp::sourceCpp("main.cpp")

tab=Fun("/home/llaniewski/Downloads/Sample 0.stl")

save(tab, file="tab.Rdata")
print(dim(tab))
print(tab[1:10,1:10])
plot(tab[,10])
plot(tab[10,])
