# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)
library(sp)           ## Data management

# Load the data
data(meuse)
d <- cbind(meuse$x, meuse$y)
y <- meuse$elev

# Find anchorpoints
a <- find_anchorpoints.lsm(d,12)
# Build the empiric variogram
vario <- variogram.lsm(y,d,a$anchorpoints,370,8,8,"gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario, "exponential", c(200,200,0.01,100))
solu
# Plot of the solutions
x11()
mypoints<-plot.lsm(model = solu, a = a, y = y, d = d, n_points = 10, points_arrangement = "random")

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d, y, d)
max(previsions$zpredicted - y)

# Test the performace of our model via cross-validation
cv.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))


y<-sic.all$data

V1<-as.numeric(sic.all$coords[,1])
V2<-as.numeric(sic.all$coords[,2])
d <- cbind(V1, V2)


# Find anchorpoints
a <- find_anchorpoints.lsm(d,8)
# Build the empiric variogram
vario <- variogram.lsm(y,d,a$anchorpoints,46*2,6,15,"gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario, "exponential", c(20,100,0.51,100),upper.bound = c(Inf,Inf,pi/2,250))
solu
plotvario(6,15,vario$empiricvariogram,34,46*2)
# Plot of the solutions
x11()
mypoints<-plot.lsm(model = solu, a = a, y = y, d = d, n_points = 10, points_arrangement = "random")

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d, y, d)
max(previsions$zpredicted - y)





