plot.lsm<-function(model, a, y, d, n_points = 4)
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
    for (j in 1:n_points)
    {
      newpoints<-rbind(newpoints, c(g$X[i]+a$width*cos(j*2*pi/n_points), g$Y[i]+a$width*sin(j*2*pi/n_points)))
    }
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
  p <- ggplot2::ggplot(dd, ggplot2::aes(x=V1, y=V2, size=y)) + ggplot2::geom_point() + ggplot2::labs(x="X", y="Y")
  p <- p + ggplot2::labs(title = "Bubble plot of the initial data", fontface = 'bold') + ggplot2::theme_light()
  print(p)
  
  
  ###ELLISSI
  ellissi<-g
  ellissi$lambda2 <- ellissi$lambda2/(ellissi$lambda1/a$width)
  ellissi$lambda1 <- a$width
  p1 <- ggplot2::ggplot(ellissi, ggplot2::aes(x=X, y=Y)) + ggforce::geom_ellipse(ggplot2::aes(x0 = X, y0 = Y, a = lambda1, b = lambda2, angle = phi), data = ellissi) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p2 <- ggplot2::ggplot(ellissi, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_segment(ggplot2::aes(x=X, y=Y, xend=X+lambda1*cos(phi), yend=Y+lambda1*sin(phi)), arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")), data = ellissi)
  p2 <- p2 + ggplot2::geom_segment(ggplot2::aes(x=X, y=Y, xend=X-lambda1*cos(phi), yend=Y-lambda1*sin(phi)), arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")), data = ellissi) + ggplot2::coord_fixed() + ggplot2::theme_light()
  print(cowplot::plot_grid(p1, p2))
  
  ###PARAMETERS
  p1 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=lambda1)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p2 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=lambda2)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p3 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=phi)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p4 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=sigma)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  title <- cowplot::ggdraw() + cowplot::draw_label("Parameters", fontface='bold')
  p <- cowplot::plot_grid(p1,p2,p3,p4)
  print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  
  
  ###FUNCTION VALUES
  predictedvalues<-predikt(y,d,model$anchorpoints,model$epsilon,model$delta,model$solutions,as.matrix(allpoints)[,1:2])
  means <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=predictedvalues$predictedmean)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
  ys <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=predictedvalues$ypredicted)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
  means<-means+ggplot2::labs(color="mean") + ggplot2::theme_light()
  ys<-ys+ggplot2::labs(color="f(*)") + ggplot2::theme_light()
  title <- cowplot::ggdraw() + cowplot::draw_label("Predicted mean and f(*)", fontface='bold')
  p <- cowplot::plot_grid(means, ys)
  print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  
  return(allpoints)
}


plot.parameters<-function(allpoints)
{
  ###PARAMETERS
  for (i in 3:dim(allpoints)[2])
  {
    par(ask=TRUE)
    p <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=allpoints[,i])) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light() + ggplot2::labs(color = colnames(allpoints)[i])
    print(p)
  }
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