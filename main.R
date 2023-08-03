
library(Rcpp)

Rcpp::sourceCpp("main.cpp")

tab = read_stl_to_tab("/home/llaniewski/Downloads/Sample 0.stl")

print(dim(tab))
print(tab[1:10,1:10])
print(attributes(tab))

# save(tab, file="tab.Rdata")

plot(tab[,10])
plot(tab[10,])

library(reticulate)
use_python("/home/llaniewski/mambaforge/bin/python3")

vtk = import("vtk")
numpy_support = import("vtk.util.numpy_support")

vtk_data = numpy_support$numpy_to_vtk(num_array=as.vector(tab))

img = vtk$vtkImageData()
img$GetPointData()$SetScalars(vtk_data)
img$SetDimensions(dim(tab)[1], dim(tab)[2], 1L)
spacing = attr(tab,"pixel")
spacing = c(spacing, mean(spacing))
img$SetSpacing(spacing[1],spacing[2],spacing[3])

writer = vtk$vtkXMLImageDataWriter()
writer$SetFileName("test.vti")
writer$SetInputData(img)
writer$Update()

