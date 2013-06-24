vcgMeshres <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    tmp <- .Call("Rmeshres",vb,it)
    return(tmp)
  }
