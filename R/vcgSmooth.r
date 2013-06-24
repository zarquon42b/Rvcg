vcgSmooth <- function(mesh,type=c("taubin","laplace","HClaplace","fujiLaplace","angWeight"),iteration=10,lambda=0.5,mu=-0.53,delta=0.1)
  {
    type <- substring(type[1],1L,1L)
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
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
    

    tmp <- .Call("Rsmooth",vb,it,iter,method,lambda,mu,delta)
    mesh$vb[1:3,] <- tmp$vb
    mesh$normals <- rbind(tmp$normals, 1)
    mesh$it <- tmp$it
    invisible(mesh)
  }

