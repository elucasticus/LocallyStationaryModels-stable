library(LocallyStationaryModels)
library(ggforce)
library(cowplot)


a<-find_anchorpoints(d,30)

r<-rawmodel(y,d,a$anchorpoints,c(200,200,0.01,100),300,8,8,"gaussian","esponenziale")
mmypoints <- plot.lsm(r,a,y,d)
