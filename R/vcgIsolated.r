vcgIsolated <- function(mesh,facenum=NULL,diameter=NULL)
  {

    if (is.null(facenum))
      facenum <- 0
    if (!is.null(diameter))
     facenum <- -1

    if (is.null(diameter))
        diameter <- 0
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    tmp <- .Call("Risolated", vb, it, diameter, facenum)
    outmesh <- list()
    class(outmesh) <- "mesh3d"
    outmesh$vb <- rbind(tmp$vb,1)
    outmesh$it <- tmp$it
    outmesh$normals <- rbind(tmp$normals, 1)
    invisible(outmesh)
  }
