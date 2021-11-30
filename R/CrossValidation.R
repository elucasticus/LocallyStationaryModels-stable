#' Cross-Validation MSE
#' 
#' @description calculate the mean squred error via cross-validation
#' @param z the vector contatining f(d)
#' @param d the matrix contatining the coordinates in which we know the value of z
#' @param anchorpoints a matrix with the coordinates of the anchorpoints which can be obtained calling find_anchorpoints.lsm
#' @param epsilon the value of epsilon regulating the kernel
#' @param n_angles the number of angles for the grid
#' @param n_intervals the number of intervals for the grid
#' @param kernel_id the type of kernel to be used
#' @param id the type of variogram to be used
#' @param initial.position the starting position to be given to the optimizer
cv.lsm <- function(z,d,anchorpoints,epsilon,n_angles,n_intervals,kernel_id, id, initial.position,n_threads = -1){
  # set the MSE to 0
  MSE=0
  
  # create the progress bar
  pb <- txtProgressBar(min = 1,         # Minimum value of the progress bar
                       max = length(z), # Maximum value of the progress bar
                       style = 3,       # Progress bar style
                       width = 50,      # Progress bar width
                       char = "=")      # Character used to create the bar
  
  for(i in 1:length(z)){
    # create a new couple (znew, dnew) deleting the i-th element from (z, d)
    znew <- z[-i]
    dnew <- d[-i,]
    
    # predict the value of f(d[i, ]) and update the MSE
    vario <- variogram.lsm(znew,dnew,anchorpoints,epsilon,n_angles,n_intervals,kernel_id,FALSE, n_threads = n_threads)
    solu <- findsolutions.lsm(vario, id, initial.position,print=FALSE, n_threads = n_threads)
    previsions <- predict.lsm(solu, rbind(d[i,]), znew, dnew,FALSE,FALSE, n_threads = n_threads)
    MSE <- MSE + (previsions$zpredicted - z[i])^2
    
    # update the progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  MSE=MSE/length(z)
  return(MSE)
}



