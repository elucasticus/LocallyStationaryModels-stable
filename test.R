# Clean the environment
rm(list = ls())

# Install and compile the package
install.packages("devtools")
library(devtools)
devtools::install_github("giacomodecarlo/LocallyStationaryModels")


# Load the libraries
library(LocallyStationaryModels)

install.packages("sp")
library(sp)           ## Only for meuse dataset

# Load the data
data(meuse)
d<- cbind(meuse$x, meuse$y)
y<-meuse$cadmium

head(d)
head(y)

# Find anchorpoints
a<-find_anchorpoints.lsm(dataset = d,n = 12)
# Build the empiric variogram
vario<-variogramlsm(y = y,d = d,anchorpoints = a$anchorpoints,epsilon = 380,n_angles = 8,n_intervals = 8,kernel_id = "gaussian")


# Plot the segments of pairs of coordinates belonging to grid cell of given index
plotgrid(d = d,grid = vario$grid,index = 3)

#Plot the anisotropic empiric variogram
plotvario(n_angles = 8,n_intervals = 8,vario$empiricvariogram,pos=22,epsilon = vario$epsilon)
#This function plots the empiric variogram for the anchor points n. "pos", along all the n_angles directions of anisotropy


# Find the parameters that minimize the nonlinear wls problem (anisotropy, variance and others)
solu<-findsolutions.lsm(vario = vario, id = "esponenziale", initial.position = c(200,200,0.01,100))
solugaussian<-findsolutions.lsm(vario = vario, id = "gaussian", initial.position = c(200,200,0.01,100,10)) #lambda1,lambda2,phi,sigma,nu
solumatern<-findsolutions.lsm(vario = vario, id = "matern", initial.position = c(200,200,0.01,100,10)) #lambda1,lambda2,phi,sigma,nu

#Variograms implemented: "exponential" "matern" "gaussian"

# Plot of the solutions
x11(height = 600, width = 800, ypos = -100, xpos = -100)
mypoints<-plot.lsm(model = solu,a = a,y = y,d = d)
mypoints<-plot.lsm(model = solugaussian,a = a,y = y,d = d)
mypoints<-plot.lsm(model = solumatern,a = a,y = y,d = d)


# Kriging on the original data
x11(height = 600, width = 800, ypos = -100, xpos = -100)
previsions <- predict.lsm(sol = solu,newpos =  d, y = y,d =  d)#solugaussian #solumatern
max(abs(previsions$zpredicted - y))
#This function can be used to perform kriging in any new set of coordinates "newpos"

