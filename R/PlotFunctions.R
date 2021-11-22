#' @brief                      generate various plots in order to better visualize the model built 
#' @param model                an object returned by findsolutions.lsm
#' @param a                    the object returned by findanchorpoints.lsm used to generate the model
#' @param z                    the vector with the values of z used to generate the model
#' @param d                    a matrix with the coordinates of the points in the original dataset used to build the model
#' @param n_points             a parameter proportional to the number of points generated to visualize the model
#' @param seed                 if points_arrangement is set to 'random', the seed used to generate the random points around each anchorpoints
#' @param points_arrangement   the arrangement of the points around each anchorpoints
plot.lsm<-function(model, a, z, d, n_points = 10, seed = 69, points_arrangement = "random", n_threads = -1, bool = TRUE)
{
  # set the seed
  set.seed(seed = seed)
  # associate each anchorpoints with the value of the parameters lambda1, lambda2, phi and sigma in its position
  aa <- as.data.frame(a$anchorpoints)
  colnames(aa) <- c("X","Y")
  s <- model$solutions
  s <- as.data.frame(s)
  colnames(s) <- c("lambda1", "lambda2", "phi", "sigma")
  g <- cbind(aa,s)
  
  # create a new dataframe to better visualize the model in the space
  # if points_arrangement = "random", generate n_points equally spaced in the angular domain at random distance (< a$width) from each anchorpoint
  # if points_arrangement = "straight", generate n_points*n_points equally spaced on straight lines around each anchorpoint
  newpoints <- data.frame(X = double(), Y = double())
  if (points_arrangement == "random")
  {
    for (i in 1:dim(g)[1])
    {
      for (j in 1:n_points)
      {
        radius <- a$width*runif(1, min=0, max=1)
        newpoints <- rbind(newpoints, c(g$X[i]+radius*cos(j*2*pi/n_points), g$Y[i]+radius*sin(j*2*pi/n_points)))
      }
    }
  }
  else if (points_arrangement == "straight")
  {
    for (i in 1:n_points)
    {
      vstep <- a$height/(n_points)
      for (j in 1:n_points)
      {
        hstep <- a$width/(n_points)
        newpoints <- rbind(newpoints, cbind(g$X-a$width+hstep*j-hstep/2, g$Y+a$height-vstep*i-vstep/2))
      }
    }
  }
  else
  {
    stop("points_arrangement is not a valid arrangement of points")
  }
  
  colnames(newpoints)<-c("X","Y")
  
  parameters<-smoothing(model$solutions,a$anchorpoints,model$delta,as.matrix(newpoints),model$kernel_id,n_threads)
  parameters<-as.data.frame(parameters)
  colnames(parameters)<-c("lambda1", "lambda2", "phi", "sigma")
  
  allpoints<-cbind(newpoints, parameters)
  
  # bubble plot of the initial data
  dd <- as.data.frame(d)
  colnames(dd) <- c("X", "Y")
  
  par(ask=TRUE)
  p <- ggplot2::ggplot(dd, ggplot2::aes(x=X, y=Y, size=z)) + ggplot2::geom_point() + ggplot2::labs(x="X", y="Y")
  p <- p + ggplot2::labs(title = "Bubble plot of the initial data", fontface = 'bold') + ggplot2::theme_light() + ggplot2::coord_fixed()
  print(p)
  
  
  # ellipses and phi
  ellissi<-g
  if (max(ellissi$lambda1) > max(ellissi$lambda2))
  {
    ellissi$lambda2 <- ellissi$lambda2/(ellissi$lambda1/a$width)
    ellissi$lambda1 <- a$width
  }
  else
  {
    ellissi$lambda1 <- ellissi$lambda1/(ellissi$lambda2/a$height)
    ellissi$lambda2 <- a$height
  }
  p1 <- ggplot2::ggplot(ellissi, ggplot2::aes(x=X, y=Y)) + ggforce::geom_ellipse(ggplot2::aes(x0 = X, y0 = Y, a = lambda1, b = lambda2, angle = phi), data = ellissi) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p2 <- ggplot2::ggplot(ellissi, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_segment(ggplot2::aes(x=X, y=Y, xend=X+a$width*cos(phi), yend=Y+a$width*sin(phi)), arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")), data = ellissi)
  p2 <- p2 + ggplot2::geom_segment(ggplot2::aes(x=X, y=Y, xend=X-a$width*cos(phi), yend=Y-a$width*sin(phi)), arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm")), data = ellissi) + ggplot2::coord_fixed() + ggplot2::theme_light()
  print(cowplot::plot_grid(p1, p2))
  
  # parameters
  p1 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=lambda1)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p2 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=lambda2)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p3 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=phi)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  p4 <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=sigma)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light()
  title <- cowplot::ggdraw() + cowplot::draw_label("Parameters", fontface='bold')
  p <- cowplot::plot_grid(p1,p2,p3,p4)
  print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  
  if(bool)
  {
    # predict and plot the mean and punctual value of z for each newpoint
    predictedvalues<-predikt(z,d,model$anchorpoints,model$epsilon,model$delta,model$solutions,as.matrix(allpoints)[,1:2],model$id,model$kernel_id,FALSE,n_threads)
    means <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=predictedvalues$predictedmean)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
    ys <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=predictedvalues$zpredicted)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
    means<-means+ggplot2::labs(color="mean") + ggplot2::theme_light()
    ys<-ys+ggplot2::labs(color="z") + ggplot2::theme_light()
    title <- cowplot::ggdraw() + cowplot::draw_label("Predicted mean and z", fontface='bold')
    p <- cowplot::plot_grid(means, ys)
    print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  }
  return(allpoints)
}

#' @brief               sequential plot of all the paramters obtained by smoothing
#' @param allpoints     a dataframe containing the coordinates and the value of the parameters
plot.parameters<-function(allpoints)
{
  for (i in 3:dim(allpoints)[2])
  {
    par(ask=TRUE)
    p <- ggplot2::ggplot(allpoints, ggplot2::aes(x=X, y=Y, color=allpoints[,i])) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed() + ggplot2::theme_light() + ggplot2::labs(color = colnames(allpoints)[i])
    print(p)
  }
}

# DA COMMENTARE
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

# DA COMMENTARE
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