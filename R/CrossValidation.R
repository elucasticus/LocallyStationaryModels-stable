cv.lsm <- function(y,d,anchorpoints,epsilon,n_angles,n_intervals,kernel_id, id, initial.position){
  MSE=0
  
  pb <- txtProgressBar(min = 1,      # Minimum value of the progress bar
                       max = length(y), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  for(i in 1:length(y)){
    
    ynew <- y[-i]
    dnew <- d[-i,]
    
    vario <-variogram.lsm(ynew,dnew,anchorpoints,epsilon,n_angles,n_intervals,kernel_id,FALSE)
    solu <-findsolutions.lsm(vario, id, initial.position,print=FALSE)
    previsions <- predict.lsm(solu, rbind(d[i,]), ynew, dnew,FALSE,FALSE)
    MSE <- MSE + (previsions$ypredicted - y[i])^2
    
    setTxtProgressBar(pb, i)
    
    
  }
  close(pb)
  MSE=MSE/length(y)
  return(MSE)
  
}



