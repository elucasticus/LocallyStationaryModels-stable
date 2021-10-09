plot.lsm<-function(model, y, d)
{
  aa<-as.data.frame(a$anchorpoints)
  colnames(aa)<-c("latitude","longitude")
  s<-model$solutions
  s<-as.data.frame(s)
  colnames(s)<-c("lambda1", "lambda2", "phi", "sigma")
  g<-cbind(aa,s)
  
  ###ELLISSI
  windows()
  p <- ggplot(aa, aes(x=latitude, y=longitude)) + geom_point(colour="black") + geom_ellipse(aes(x0 = latitude, y0 = longitude, a = lambda1, b = lambda2, angle = phi), data = g) + coord_fixed()
  print(p)
  
  ###DATI INIZIALI
  dd <- as.data.frame(d)
  windows()
  p <- ggplot(dd, aes(x=V1, y=V2, color=y)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  print(p)
  
  windows()
  p <- ggplot(dd, aes(x=V1, y=V2, size=y)) + geom_point()
  print(p)
}