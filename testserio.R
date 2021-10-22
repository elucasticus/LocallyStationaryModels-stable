rm(list = ls())

library(LocallyStationaryModels)
library(ggforce)
library(cowplot)
library(sp)           ## Data management
data(meuse)
d <- cbind(meuse$x, meuse$y)
y<-meuse$cadmium

dnew=(d[!(d[,1]>180000&d[,2]<330500),])
ynew=y[!(d[,1]>180000&d[,2]<330500)]

a<-find_anchorpoints(dnew,30)
vario<-variogramlsm(ynew,dnew,a$anchorpoints,290,8,8,"gaussian")
solu<-findsolutions.lsm(vario, "esponenziale", c(200,200,0.01,100))

x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints<-plot.lsm(solu,a,ynew,dnew)

x11(height = 600, width = 800, ypos = -100, xpos = -100)
previsions <- predict.lsm(solu, dnew, ynew, dnew)
max(previsions$ypredicted - ynew)
