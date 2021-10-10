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
    newpoints<-rbind(newpoints, c(g$latitude[i]+a$width/2, g$longitude[i]))
    newpoints<-rbind(newpoints, c(g$latitude[i]-a$width/2, g$longitude[i]))
    newpoints<-rbind(newpoints, c(g$latitude[i], g$longitude[i]+a$height/2))
    newpoints<-rbind(newpoints, c(g$latitude[i], g$longitude[i]-a$height/2))
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
  windows()
  p <- ggplot(dd, aes(x=V1, y=V2, color=y)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  print(p)
  
  windows()
  p <- ggplot(dd, aes(x=V1, y=V2, size=y)) + geom_point()
  print(p)
  
  
  ###ELLISSI
  windows()
  p <- ggplot(g, aes(x=latitude, y=longitude)) + geom_point(colour="black") + geom_ellipse(aes(x0 = latitude, y0 = longitude, a = lambda1, b = lambda2, angle = phi), data = g) + coord_fixed()
  print(p)
  
  ###PARAMETERS
  windows()
  p <- ggplot(allpoints, aes(x=latitude, y=longitude, color=sigma)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  print(p)
}