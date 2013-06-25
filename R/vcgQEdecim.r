vcgQEdecim <- function(mesh,tarface=NULL,percent=NULL,edgeLength=NULL, topo=TRUE,quality=TRUE,bound=TRUE, optiplace = TRUE, scaleindi = TRUE, normcheck = FALSE, safeheap =FALSE, qthresh=0.1, boundweight = 0.5, normalthr = pi/2)
  {
    doit <- TRUE
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    dimit <- ncol(it)
    outmesh <- list()
    class(outmesh) <- "mesh3d"
    
    if (is.null(tarface) && is.null(percent)&& is.null(edgeLength))
      stop("please enter decimation option")
    if (!is.null(percent))
      {
        if (percent <= 0 || percent > 1)
          stop ("percent must be between 0 and 1")
        tarface <- floor(percent*dimit)
      }
    if (!is.null(edgeLength))
      {
        currEL <- vcgMeshres(mesh)$res
        if (currEL >= edgeLength)
          warning("edges already shorter than required - nothing to do")
        coef <- (currEL/edgeLength)^2
        tarface <- floor(coef*dimit)
        
      }
    ##concatenate parameters
    boolparams <- c( topo, quality, bound, optiplace, scaleindi, normcheck, safeheap)
    doubleparams <- c(qthresh, boundweight, normalthr)
###tmp <- .C("RQEdecim",vb,ncol(vb),it,ncol(it),tarface,vb,as.integer(topo),as.integer(quality),as.integer(bound))
    tmp <- .Call("RQEdecim", vb, it, tarface, boolparams, doubleparams)
    outmesh$vb <- rbind(tmp$vb,1)
    
    outmesh$it <- tmp$it
    outmesh$normals <- rbind(tmp$normals, 1)
    #outmesh <- adnormals(outmesh)
    if(!is.null(edgeLength))
    cat(paste("Mean Edge length is",vcgMeshres(outmesh)$res,"\n"))
    return(outmesh)
  }
