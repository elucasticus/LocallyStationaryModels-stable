findsolutions.lsm<-function(vario, id, initial.position)
{
  return(findsolutionslsm(vario$anchorpoints,vario$empiricvariogram,vario$squaredweigths,vario$mean.x, vario$mean.y, id, initial.position,vario$epsilon))
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