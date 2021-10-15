library(LocallyStationaryModels)
library(ggforce)
library(cowplot)
library(sp)           ## Data management
data(meuse)
y<-meuse$cadmium


a<-find_anchorpoints(d,40)

r<-rawmodel(y,d,a$anchorpoints,c(200,200,0.01,100),350,8,8,"gaussian","esponenziale")
x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints <- plot.lsm(r,a,y,d)

griglia2<-buildgrid(y,d,a$anchorpoints,500,8,8,"gaussian")
plotgrid(d,griglia$grid ,9)

plotvario(8,8,griglia$empiricvariogram,5,500)


###IN ALTERNATIVA
vario<-variogramlsm(y,d,a$anchorpoints,350,8,8,"gaussian")
solu<-findsolutionslsm(vario$anchorpoints,vario$empiricvariogram,vario$squaredweigths,vario$mean.x, vario$mean.y, "esponenziale", c(200,200,0.01,100),vario$epsilon)
mypoints<-plot.lsm(solu,a,y,d)

###IN ALTERNATIVA (DA SISTEMARE)
vario<-variogramlsm(y,d,a$anchorpoints,350,8,8,"gaussian")
solu<-findsolutions.lsm(vario, "esponenziale", c(200,200,0.01,100))
previsions <- predict.lsm(solu, d, y, d, TRUE)
