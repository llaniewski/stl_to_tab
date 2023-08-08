
library(Rcpp)
library(reticulate)
use_python("/home/llaniewski/mambaforge/bin/python3")

vtk = import("vtk")
numpy_support = import("vtk.util.numpy_support")

Rcpp::sourceCpp("main.cpp")

#samples=c("sample_0","sample_1_mag10","sample_1_mag40")
samples=c("sample_2_mag10","sample_2_mag40")

for (sample in samples) {

print(sample)

tab = read_stl_to_tab(paste0("samples/",sample,".stl"))

print(dim(tab))
print(tab[1:10,1:10])
print(attributes(tab))

save(tab, file=paste0(sample,".Rdata"))

plot(tab[,10])
plot(tab[10,])

vtk_data = numpy_support$numpy_to_vtk(num_array=as.vector(tab))

img = vtk$vtkImageData()
img$GetPointData()$SetScalars(vtk_data)
img$SetDimensions(dim(tab)[1], dim(tab)[2], 1L)
spacing = attr(tab,"pixel")
spacing = c(spacing, mean(spacing))
img$SetSpacing(spacing[1],spacing[2],spacing[3])

writer = vtk$vtkXMLImageDataWriter()
writer$SetFileName(paste0(sample,".vti"))
writer$SetInputData(img)
writer$Update()

}


