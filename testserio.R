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
vario <- variogram.lsm(y,d,a$anchorpoints,350,8,8,"gaussian")
# Find the solutions
solu <- findsolutions.lsm(vario, "matern", c(200,200,0.01,100))

# Plot of the solutions
x11()
mypoints<-plot.lsm(model = solu, a = a, y = y, d = d, n_points = 10, points_arrangement = "random")

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d, y, d)
max(previsions$ypredicted - y)

# Test the performace of our model via cross-validation
cv.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))
