library(LocallyStationaryModels)
library(ggforce)
library(cowplot)
library(sp)           ## Data management
data(meuse)
k<-meuse$cadmium


a<-find_anchorpoints(d,40)

r<-rawmodel(k,d,a$anchorpoints,c(200,200,0.01,100),350,8,8,"gaussian","esponenziale")
mmypoints <- plot.lsm(r,a,k,d)





griglia2<-buildgrid(y,d,a2$anchorpoints,500,8,8,"gaussian")
plotgrid(d,griglia$grid ,9)

plotvario(8,8,griglia$empiricvariogram,5,500)


