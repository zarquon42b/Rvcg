vcgIsolated <- function(mesh,facenum=NULL,diameter=NULL)
  {

    if (is.null(facenum))
      facenum <- 0
    if (!is.null(diameter))
     facenum <- -1
    outmesh <- list()
    class(outmesh) <- "mesh3d"
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    storage.mode(it) <- "integer"
    storage.mode(facenum) <- "integer"
    storage.mode(diameter) <- "double"
    tmp <- .C("Risolated",vb,ncol(vb),it,ncol(it),diameter,vb,facenum)
   
    outmesh$vb <- rbind(tmp[[1]][,1:tmp[[2]]],1)
    outmesh$it <- tmp[[3]][,1:(tmp[[4]])]+1
    outmesh$normals <- rbind(tmp[[6]][,1:(tmp[[2]])], 1)
    invisible(outmesh)
  }
