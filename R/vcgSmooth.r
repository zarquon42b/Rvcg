vcgSmooth <- function(mesh,type=c("taubin","laplace","HClaplace"),iteration=10,sign=F)
  {
    type <- substring(type[1],1L,1L)
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    method <- 0
    if (type == "l")
      {
        method <- 1
      }
     else if (type == "H" || type == "h")
      {
        method <- 2
      }
    iter=iteration
    storage.mode(it) <- "integer"
    storage.mode(method) <- "integer"
    storage.mode(iter) <- "integer"
    normals <- vb

    tmp <- .C("Rsmooth",vb,ncol(vb),it,ncol(it),iter,method,normals)
    mesh$vb[1:3,] <- tmp[[1]]
    mesh$normals <- rbind(tmp[[7]],1)
    invisible(mesh)
  }

