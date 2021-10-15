findsolutions.lsm<-function(vario, id, initial.position)
{
  return(findsolutionslsm(vario$anchorpoints,vario$empiricvariogram,vario$squaredweigths,vario$mean.x, vario$mean.y, id, initial.position,vario$epsilon))
}