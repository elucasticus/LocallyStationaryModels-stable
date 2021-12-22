# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)

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
# Plot of the solutions
x11()
mypoints<-plot.lsm(model = solu, a = a, z = y, d = d, n_points = 3, points_arrangement = "straight", kriging = FALSE, 
                   ellipse_scale = 2, arrow_scale = 1.5)

# Kriging on the original data
x11()
previsions <- predict.lsm(solu, d)
max(previsions$zpredicted - y)

# Test the performace of our model via cross-validation
cv.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))
