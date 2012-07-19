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
    storage.mode(clost) <- "double"
    sign <- as.integer(sign)
    tmp <- .C("Rclost",vb,ncol(vb),it,ncol(it),clost,clostDim,normals,dis,sign)
    x$vb[1:3,] <- tmp[[5]]
    x$normals <- rbind(tmp[[7]],1)
    chcknorm <- which(is.nan(x$normals))
    if (length(chcknorm) > 0)
      x$normals[chcknorm] <- 0
                      
    x$quality <- tmp[[8]]

    invisible(x)
  }
    
