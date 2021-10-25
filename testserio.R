# Clean the environment
rm(list = ls())

# Load the libraries
library(LocallyStationaryModels)
library(ggforce)
library(cowplot)
library(sp)           ## Data management

# Load the data
data(meuse)
d <- cbind(meuse$x, meuse$y)
y<-meuse$cadmium

# Delete problematic data points
dnew=(d[!(d[,1]>180000&d[,2]<330500),])
ynew=y[!(d[,1]>180000&d[,2]<330500)]


head(dnew)
head(ynew)

# Find anchorpoints
a<-find_anchorpoints(dnew,30)
# Build the empiric variogram
vario<-variogramlsm(ynew,dnew,a$anchorpoints,290,8,8,"gaussian")
# Find the solutions
solu<-findsolutions.lsm(vario, "esponenziale", c(200,200,0.01,100))

# Plot of the solutions
x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints<-plot.lsm(solu,a,ynew,dnew)

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
