plot.lsm<-function(model, a, y, d)
{
  aa<-as.data.frame(a$anchorpoints)
  colnames(aa)<-c("X","Y")
  s<-model$solutions
  s<-as.data.frame(s)
  colnames(s)<-c("lambda1", "lambda2", "phi", "sigma")
  g<-cbind(aa,s)
  
  newpoints <- data.frame(X = double(), Y = double())
  for (i in 1:dim(g)[1])
  {
    #newpoints<-rbind(newpoints, c(g$X[i]+a$width/2, g$Y[i]))
    #newpoints<-rbind(newpoints, c(g$X[i]-a$width/2, g$Y[i]))
    #newpoints<-rbind(newpoints, c(g$X[i], g$Y[i]+a$height/2))
    #newpoints<-rbind(newpoints, c(g$X[i], g$Y[i]-a$height/2))
    newpoints<-rbind(newpoints, c(g$X[i]+a$width/sqrt(2), g$Y[i]+a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$X[i]+a$width/sqrt(2), g$Y[i]-a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$X[i]-a$width/sqrt(2), g$Y[i]+a$height/sqrt(2)))
    newpoints<-rbind(newpoints, c(g$X[i]-a$width/sqrt(2), g$Y[i]-a$height/sqrt(2)))
  }
  colnames(newpoints)<-c("X","Y")
  
  parameters<-smoothing(model$solutions,a$anchorpoints,model$delta,as.matrix(newpoints))
  parameters<-as.data.frame(parameters)
  colnames(parameters)<-c("lambda1", "lambda2", "phi", "sigma")
  
  allpoints<-cbind(newpoints, parameters)
  #allpoints<-rbind(g, newpos)
  
  ###DATI INIZIALI
  dd <- as.data.frame(d)
  #windows()
  #p <- ggplot(dd, aes(x=V1, y=V2, color=y)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  #print(p)
  
  par(ask=TRUE)
  p <- ggplot(dd, aes(x=V1, y=V2, size=y)) + geom_point() + labs(x="X", y="Y")
  p <- p + labs(title = "Bubble plot of the initial data", fontface = 'bold') + theme_light()
  print(p)
  
  
  ###ELLISSI
  ellissi<-g
  ellissi$lambda2 <- ellissi$lambda2/(ellissi$lambda1/a$width)
  ellissi$lambda1 <- a$width
  p <- ggplot(ellissi, aes(x=X, y=Y)) + geom_ellipse(aes(x0 = X, y0 = Y, a = lambda1, b = lambda2, angle = phi), data = ellissi) + coord_fixed() + theme_light()
  print(p)
  
  ###PARAMETERS
  p1 <- ggplot(allpoints, aes(x=X, y=Y, color=lambda1)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed() + theme_light()
  p2 <- ggplot(allpoints, aes(x=X, y=Y, color=lambda2)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed() + theme_light()
  p3 <- ggplot(allpoints, aes(x=X, y=Y, color=phi)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed() + theme_light()
  p4 <- ggplot(allpoints, aes(x=X, y=Y, color=sigma)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed() + theme_light()
  title <- ggdraw() + draw_label("Parameters", fontface='bold')
  p <- plot_grid(p1,p2,p3,p4)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  
  
  ###FUNCTION VALUES
  predictedvalues<-predikt(y,d,model$anchorpoints,model$epsilon,model$delta,model$solutions,as.matrix(allpoints)[,1:2])
  means <- ggplot(allpoints, aes(x=X, y=Y, color=predictedvalues$predictedmean)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  ys <- ggplot(allpoints, aes(x=X, y=Y, color=predictedvalues$ypredicted)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
  means<-means+labs(color="mean") + theme_light()
  ys<-ys+labs(color="f(*)") + theme_light()
  title <- ggdraw() + draw_label("Predicted mean and f(*)", fontface='bold')
  p <- plot_grid(means, ys)
  print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  
  return(allpoints)
}



plotgrid<-function(d,grid,index){
  
  plot(d,xlab="Latitude",ylab="Longitude")
  grid=as.data.frame(grid)
  
  for (i in 1:dim(d)[1]){
    for (j in i:dim(d)[1]){
      
      if(grid[i,j]==index){
        segments(d[i,1],d[i,2],d[j,1],d[j,2])
                 
      }
      
    }
  }
  
}


plotvario<-function(n_angles,n_intervals,empvariogram,pos,epsilon){
  b=2*epsilon
  diminterval = b/n_intervals
  
  coordnormh = numeric(n_intervals)
  for ( i in 1:n_intervals){
    coordnormh[i] = diminterval/2 + (i-1)*diminterval
  }
  
  par(ask=TRUE)
  for (i in 1:n_angles){
    plot(coordnormh , empvariogram[(n_intervals*(i-1)+1):(n_intervals*i),pos] , xlab=paste("Distance",
      "  (",
      as.character(180*(i-1)/n_angles),"\u00B0 -",
      as.character(180*i/n_angles),"\u00B0 )"),
         ylab= paste("Empiric Anisotropic Variogram"))
  }
  
  
  
}