plot.lsm<-function(model, a, y, d)
{
  aa<-as.data.frame(a$anchorpoints)
  colnames(aa)<-c("latitude","longitude")
  s<-model$solutions
  s<-as.data.frame(s)
  colnames(s)<-c("lambda1", "lambda2", "phi", "sigma")
  g<-cbind(aa,s)
  
  newpoints <- data.frame(latitude = double(), longitude = double())
  for (i in 1:dim(g)[1])
  {
    #newpoints<-rbind(newpoints, c(g$latitude[i]+a$width/2, g$longitude[i]))
    #newpoints<-rbind(newpoints, c(g$latitude[i]-a$width/2, g$longitude[i]))
    #newpoints<-rbind(newpoints, c(g$latitude[i], g$longitude[i]+a$height/2))
    #newpoints<-rbind(newpoints, c(g$latitude[i], g$longitude[i]-a$height/2))
    newpoints<-rbind(newpoints, c(g$latitude[i]+a$width/sqrt(2), g$longitude[i]+a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$latitude[i]+a$width/sqrt(2), g$longitude[i]-a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$latitude[i]-a$width/sqrt(2), g$longitude[i]+a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$latitude[i]-a$width/sqrt(2), g$longitude[i]-a$height/sqrt(2)))
  }
  colnames(newpoints)<-c("latitude","longitude")
  
  parameters<-smoothing(model$solutions,a$anchorpoints,model$delta,as.matrix(newpoints))
  parameters<-as.data.frame(parameters)
  colnames(parameters)<-c("lambda1", "lambda2", "phi", "sigma")
  
  newpos<-cbind(newpoints, parameters)
  allpoints<-rbind(g, newpos)
  
  ###DATI INIZIALI
  dd <- as.data.frame(d)
  #windows()
  #p <- ggplot(dd, aes(x=V1, y=V2, color=y)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  #print(p)
  
  x11(ypos=-100, xpos=-100)
  p <- ggplot(dd, aes(x=V1, y=V2, size=y)) + geom_point()
  print(p)
  
  
  ###ELLISSI
  x11(ypos=-100,xpos=-100)
  p <- ggplot(g, aes(x=latitude, y=longitude)) + geom_ellipse(aes(x0 = latitude, y0 = longitude, a = lambda1, b = lambda2, angle = phi), data = g) + coord_fixed()
  print(p)
  
  ###PARAMETERS
  x11(height=600, width=600, ypos=-100,xpos=-100)
  p1 <- ggplot(allpoints, aes(x=latitude, y=longitude, color=lambda1)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  p2 <- ggplot(allpoints, aes(x=latitude, y=longitude, color=lambda2)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  p3 <- ggplot(allpoints, aes(x=latitude, y=longitude, color=phi)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  p4 <- ggplot(allpoints, aes(x=latitude, y=longitude, color=sigma)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  print(plot_grid(p1,p2,p3,p4,labels="AUTO"))
  
  
  ###FUNCTION VALUES
  predictedvalues<-predikt(y,d,model$anchorpoints,model$epsilon,model$delta,model$solutions,as.matrix(allpoints)[,1:2])
  x11(ypos=-100, xpos=-100, height=300, width=800)
  means <- ggplot(allpoints, aes(x=latitude, y=longitude, color=predictedvalues$predictedmean)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  ys <- ggplot(allpoints, aes(x=latitude, y=longitude, color=predictedvalues$ypredicted)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  print(plot_grid(means, ys,labels="AUTO"))
}