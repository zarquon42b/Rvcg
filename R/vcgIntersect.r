vcgIntersect <- function(x,mesh)
{
  
  vb <- mesh$vb[1:3,]
  it <- mesh$it - 1
  dimit <- dim(it)[2]
  dimvb <- dim(vb)[2]
  storage.mode(it) <- "integer"
  if (is.null(x$normals))
    adnormals(x)
  clost <- x$vb[1:3,]
  normals <- x$normals[1:3,]
  clostDim <- ncol(clost)
  dis <- rep(0,clostDim)
  hit <- dis
  storage.mode(hit) <- "integer"
  tmp <- .C("Rintersect",vb,ncol(vb),it,ncol(it),clost,clostDim,normals,dis,hit)
  x$vb[1:3,] <- tmp[[5]]
  x$normals <- rbind(tmp[[7]],1)
  x$quality <- tmp[[9]]
  x$distance <- tmp[[8]]
  return(x)
}
