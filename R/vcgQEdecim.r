vcgQEdecim <- function(mesh,tarface=NULL,percent=NULL)
  {
    
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    outmesh <- list()
    class(outmesh) <- "mesh3d"
    storage.mode(it) <- "integer"
    
    if (is.null(tarface) && is.null(percent))
      stop("please enter decimation option")
    if (!is.null(percent))
      {
        if (percent <= 0 || percent > 1)
          stop ("percent must be between 0 and 1")
        tarface <- floor(percent*dimit)
      }
    storage.mode(tarface) <- "integer"
    tmp <- .C("RQEdecim",vb,ncol(vb),it,ncol(it),tarface,vb)
    outmesh$vb <- rbind(tmp[[1]][,1:tmp[[2]]],1)
    
    outmesh$it <- tmp[[3]][,1:(tmp[[4]])]+1
    outmesh$normals <- rbind(tmp[[6]][,1:(tmp[[2]])], 1)
    #outmesh <- adnormals(outmesh)
    
    return(outmesh)
  }
