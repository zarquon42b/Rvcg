vcgSmooth <- function(mesh,type=c("taubin","laplace","HClaplace","fujiLaplace","angWeight"),iteration=10,lambda=0.5,mu=-0.53,delta=0.1)
  {
    type <- substring(type[1],1L,1L)
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    method <- 0
    if (type == "l" || type == "L")
      {
        method <- 1
      }
     else if (type == "H" || type == "h")
      {
        method <- 2
      }
    else if (type == "f" || type == "F")
      {
        method <- 3
      }
    else if (type == "a" || type == "A")
      {
        method <- 4
      }
    iter=iteration
    storage.mode(it) <- "integer"
    storage.mode(method) <- "integer"
    storage.mode(iter) <- "integer"
    normals <- vb

    tmp <- .C("Rsmooth",vb,ncol(vb),it,ncol(it),iter,method,normals,lambda,mu,delta)
    mesh$vb[1:3,] <- tmp[[1]]
    mesh$normals <- rbind(tmp[[7]],1)
    invisible(mesh)
  }

