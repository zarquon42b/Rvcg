vcgMeshres <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    storage.mode(it) <- "integer"
    
    tmp <- .C("Rmeshres",vb,ncol(vb),it,ncol(it),0)
    return(tmp[[5]])
  }
