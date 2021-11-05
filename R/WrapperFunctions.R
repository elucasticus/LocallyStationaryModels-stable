#' @brief                     for each anchorpoints solves a problem of nonlinear optimization and returns the results
#' @param vario               a variogram obtained using variogram.lsm()
#' @param id                  the type of variogram to be used
#' @param initial.position    the starting position to be given to the optimizer
#' @param bool                if set to TRUE removes the anchorpoints which cause troubles to the optimizer
findsolutions.lsm<-function(vario, id, initial.position, bool = FALSE, print = TRUE)
{
  if(grepl("maternNuFixed", id, fixed = TRUE))
  {
    id_check <- "maternNuFixed"
  }
  else
  {
    id_check <- id
  }
  load("variogramfunctions.Rdata")
  if(length(initial.position) != variogramfunctions$n_parameters[which(variogramfunctions$name == id_check)])
  {
    stop("wrong number of initial parameters")
  }
  result <- findsolutionslsm(vario$anchorpoints,vario$empiricvariogram,vario$squaredweigths,vario$mean.x, vario$mean.y, id, vario$kernel_id, initial.position,vario$epsilon,print)
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
  result$id <- id
  result$kernel_id <- vario$kernel_id
  class(result) <- "lsm"
  return(result)
}

#' @brief         for each couple of coordinates in newpos predict the mean and punctual value of f(*)
#' @param sol     an object of type lsm obtained by calling findsolutions.lsm
#' @param newpos  a matrix with the coordinates of the points where to evaluate f(*)
#' @param y       the vector y used to generate the solutions
#' @param d       the matrix d used to generate the solutions
#' @param bool    if set to TRUE plot the solutions
predict.lsm<-function(sol, newpos, y, d, bool = TRUE, print = TRUE)
{
  if(length(y) != dim(d)[1])
  {
    print("The length of y and the number or rows of d do not coincide")
  }
  predictedvalues <- predikt(y,d,sol$anchorpoints,sol$epsilon,sol$delta,sol$solutions,newpos,sol$id,sol$kernel_id,print)
  if (bool)
  {
    newpos <- as.data.frame(newpos)
    colnames(newpos) <- c("X", "Y")
    means <- ggplot2::ggplot(newpos, ggplot2::aes(x=X, y=Y, color=predictedvalues$predictedmean)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
    ys <- ggplot2::ggplot(newpos, ggplot2::aes(x=X, y=Y, color=predictedvalues$ypredicted)) + ggplot2::geom_point() + ggplot2::scale_color_gradientn(colours = rainbow(5)) + ggplot2::coord_fixed()
    means<-means+ggplot2::labs(color="mean") + ggplot2::theme_light()
    ys<-ys+ggplot2::labs(color="f(*)") + ggplot2::theme_light()
    title <- cowplot::ggdraw() + cowplot::draw_label("Predicted mean and f(*)", fontface='bold')
    p <- cowplot::plot_grid(means, ys)
    print(cowplot::plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)))
  }
  return(predictedvalues)
}

#' @brief          given a dataset find the corresponding equally spaced anchorpoints
#' @param dataset  a dataset full of coordinates
#' @param n        a parameter proportional to the number of anchorpoints
#' @param bool     if set to true plot the original dataset and the anchorpoints
find_anchorpoints.lsm<-function(dataset, n, bool = TRUE)
{
  result <- find_anchorpoints(dataset, n)
  if (bool)
  {
    D <- as.data.frame(dataset)
    A <- as.data.frame(result)
    colnames(D) <- c("X", "Y")
    colnames(A) <- c("X", "Y")
    p <- ggplot2::ggplot(data = D, ggplot2::aes(x=X, y=Y)) + ggplot2::geom_point() + ggplot2::geom_point(data = A, ggplot2::aes(x=X, y=Y, color="red"), shape = 3) + ggplot2::theme_light() + ggplot2::theme(legend.position = "none") + ggplot2::coord_fixed()
    print(p)
  }
  return(result)
}

#' @brief               compute the sample variogram in the anchorpoints
#' @param y             the vector contatining f(d)
#' @param d             the matrix contatining the coordinates in which we know the value of f(*)
#' @param anchorpoints  a matrix with the coordinates of the anchorpoints which can be obtained calling find_anchorpoints.lsm
#' @param epsilon       the value of epsilon regulating the kernel
#' @param n_angles      the number of angles for the grid
#' @param n_intervals   the number of intervals for the grid
#' @param kernel_id     the type of kernel to be used
variogram.lsm <- function(y, d, anchorpoints, epsilon, n_angles, n_intervals, kernel_id, print=TRUE)
{
  if(length(y) != dim(d)[1])
  {
    print("The length of y and the number or rows of d do not coincide")
  }
  vario <- variogramlsm(y, d, anchorpoints, epsilon, n_angles, n_intervals, kernel_id,print)
  vario$kernel_id <- kernel_id
  return(vario)
}
