vcgSmooth <- function(mesh,type=c("taubin","laplace"),iteration=10,sign=F)
  {
    type <- type[1]
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    method <- 0
    if (type == "laplace")
      {
        method <- 1
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
vcgClost <- function(x,mesh,sign=TRUE)
  
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it - 1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    storage.mode(it) <- "integer"

    if (is.matrix(x))
      {
        clost <- t(x)
        x <- list()
        x$vb <- clost
      }
    else
      {
        clost <- x$vb[1:3,]
      }
    normals <- clost
    clostDim <- ncol(clost)
    dis <- rep(0,clostDim)
    sign <- as.integer(sign)
    tmp <- .C("Rclost",vb,ncol(vb),it,ncol(it),clost,clostDim,normals,dis,sign)
    x$vb[1:3,] <- tmp[[5]]
    x$normals <- rbind(tmp[[7]],1)
    x$quality <- tmp[[8]]

    invisible(x)
  }
    
