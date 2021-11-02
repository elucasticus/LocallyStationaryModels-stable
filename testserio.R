# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)
library(sp)           ## Data management

# Load the data
data(meuse)
d <- cbind(meuse$x, meuse$y)
y<-meuse$cadmium



head(dnew)
head(ynew)

cv.lsm(y,d,a$anchorpoints,350,8,8,"gaussian","exponential", c(200,200,0.01,100))

# Find anchorpoints
a<-find_anchorpoints.lsm(dnew,12)
# Build the empiric variogram
vario<-variogram.lsm(ynew,dnew,a$anchorpoints,350,8,8,"gaussian")
# Find the solutions
solu<-findsolutions.lsm(vario, "exponential", c(200,200,0.01,100))

# Plot of the solutions
x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints<-plot.lsm(model = solu, a = a, y = ynew,d = dnew, n_points = 10, points_arrangement = "random")

# Kriging on the original data
x11(height = 600, width = 800, ypos = -100, xpos = -100)
previsions <- predict.lsm(solu, dnew, ynew, dnew)
max(previsions$ypredicted - ynew)



# Test the performace of our model
n <- sample(1:100, 10)
dtrain <- dnew[-n, ]
ytrain <- ynew[-n]
dvalidate <- dnew[n, ]
yvalidate <- ynew[n]

head(dtrain)
head(ytrain)
head(dvalidate)
head(yvalidate)

# Find anchorpoints
atrain <- find_anchorpoints(dtrain,30)
# Build the empiric variogram
variotrain <- variogramlsm(ytrain,dtrain,atrain$anchorpoints,290,8,8,"gaussian")
# Find the solutions
solutrain <- findsolutions.lsm(variotrain, "esponenziale", c(200,200,0.01,100))

# Plot of the solutions
x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints<-plot.lsm(solutrain,atrain,ytrain,dtrain)

# Kriging on the validate dataset
x11(height = 600, width = 800, ypos = -100, xpos = -100)
previsionstrain <- predict.lsm(solu, dvalidate, ytrain, dtrain)









# Build the empiric variogram
v1<-variogramlsm(ynew,dnew,a$anchorpoints,300,8,8,"gaussian")
# Find the solutions
s1<-findsolutions.lsm(v1, "esponenziale", c(200,200,0.01,100))
# Build the empiric variogram
v2<-variogramlsm(ynew,dnew,a$anchorpoints,3000,8,8,"gaussian")
# Find the solutions
s2<-findsolutions.lsm(v2, "esponenziale", c(200,200,0.01,100))


