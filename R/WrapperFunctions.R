findsolutions.lsm<-function(vario, id, initial.position, bool = FALSE)
{
  result <- findsolutionslsm(vario$anchorpoints,vario$empiricvariogram,vario$squaredweigths,vario$mean.x, vario$mean.y, id, initial.position,vario$epsilon)
  if (bool)
  {
    for (i in 1:dim(result$solutions)[1])
    {
      if (norm((result$solutions[i, ]-initial.position), type="2") < 1e-12)
      {
        result$solutions <- result$solutions[-i,]
        result$anchorpoints <- result$anchorpoints[-i, ]
      }
    }
  }
  class(result) <- "lsm"
  return(result)
}


predict.lsm<-function(sol, newpos, y, d, bool = TRUE)
{
  predictedvalues <- predikt(y,d,sol$anchorpoints,sol$epsilon,sol$delta,sol$solutions,newpos)
  if (bool)
  {
    newpos <- as.data.frame(newpos)
    colnames(newpos) <- c("X", "Y")
    means <- ggplot(newpos, aes(x=X, y=Y, color=predictedvalues$predictedmean)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
    ys <- ggplot(newpos, aes(x=X, y=Y, color=predictedvalues$ypredicted)) + geom_point() + scale_color_gradientn(colours = rainbow(5)) + coord_fixed()
    means<-means+labs(color="mean") + theme_light()
    ys<-ys+labs(color="f(*)") + theme_light()
    title <- ggdraw() + draw_label("Predicted mean and f(*)", fontface='bold')
    p <- plot_grid(means, ys)
    print(plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  }
  return(predictedvalues)
}


find_anchorpoints.lsm<-function(dataset, n, bool = TRUE)
{
  result <- find_anchorpoints(dataset, n)
  if (bool)
  {
    D <- as.data.frame(dataset)
    A <- as.data.frame(result)
    colnames(D) <- c("X", "Y")
    colnames(A) <- c("X", "Y")
    p <- ggplot(data = D, aes(x=X, y=Y)) + geom_point() + geom_point(data = A, aes(x=X, y=Y, color="red"), shape = 3) + labs(color="anchorpoints") + theme_light() + coord_fixed()
    print(p)
  }
  return(result)
}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Package LocallyStationaryModels version 1.0\nThis package uses OpenMP to speed up computations\nWritten by Luca Crippa and Giacomo de Carlo")
}