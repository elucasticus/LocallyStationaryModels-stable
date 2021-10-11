library(ggforce)
library(cowplot)


a<-find_anchorpoints(d,8)

r<-rawmodel(y,d,a$anchorpoints,c(200,200,0.01,100),550008,8,12,"gaussian","esponenziale")
mmypoints <- plot.lsm(r,a,y,d)
