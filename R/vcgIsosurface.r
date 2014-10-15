#' Create Isosurface from 3D-array
#'
#' Create Isosurface from 3D-array using Marching Cubes algorithm
#'
#' @param vol an integer valued 3D-array
#' @param lower numeric:lower threshold
#' @param upper numeric: upper threshold
#' @param spacing numeric 3D-vector: specifies the voxel dimensons in x,y,z direction.
#' @param origin numeric 3D-vector: origin of the original data set, will transpose the mesh onto that origin.
#' @param threshold threshold of intersecting the cube (default is 0.5).
#' @return returns a triangular mesh of class "mesh3d"
#' @examples
#' #this is the example from the package "misc3d"
#' x <- seq(-2,2,len=50)
#' g <- expand.grid(x = x, y = x, z = x)
#' v <- array(g$x^4 + g$y^4 + g$z^4, rep(length(x),3))
#' storage.mode(v) <- "integer"
#' \dontrun{
#' mesh <- vcgIsosurface(v,lower=1)
#' require(rgl)
#' wire3d(mesh)
#' ##now smooth it a little bit
#' wire3d(vcgSmooth(mesh,"HC",iteration=3),col=3)
#' }
#' @export
vcgIsosurface <- function(vol,lower=min(vol),upper=max(vol),spacing=NULL, origin=NULL,threshold=0.5) {
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    if (threshold < 0 || threshold >= 1)
        stop("threshold must be >=0 and < 1")
    lower <- as.numeric(lower)
    upper <- as.numeric(upper)
    storage.mode(vol) <- "integer"
    volmesh <- .Call("RMarchC",vol,lower,upper,threshold)
    volmesh$vb <- rbind(volmesh$vb,1)
    volmesh$it <- volmesh$it
    
    class(volmesh) <- "mesh3d"
    if (!is.null(spacing))
        volmesh$vb[1:3,] <- volmesh$vb[1:3,]*spacing
    
    if (!is.null(origin))
            volmesh$vb[1:3,] <- volmesh$vb[1:3,]+origin
    
    return(volmesh)
}

