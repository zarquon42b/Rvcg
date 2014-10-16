#' Create Isosurface from 3D-array
#'
#' Create Isosurface from 3D-array using Marching Cubes algorithm
#'
#' @param vol an integer valued 3D-array
#' @param spacing numeric 3D-vector: specifies the voxel dimensons in x,y,z direction.
#' @param origin numeric 3D-vector: origin of the original data set, will transpose the mesh onto that origin.
#' @param threshold threshold for creating the surface 
#' @return returns a triangular mesh of class "mesh3d"
#' @examples
#' #this is the example from the package "misc3d"
#' x <- seq(-2,2,len=50)
#' g <- expand.grid(x = x, y = x, z = x)
#' v <- array(g$x^4 + g$y^4 + g$z^4, rep(length(x),3))
#' storage.mode(v) <- "integer"
#' \dontrun{
#' mesh <- vcgIsosurface(v,threshold=10)
#' require(rgl)
#' wire3d(mesh)
#' ##now smooth it a little bit
#' wire3d(vcgSmooth(mesh,"HC",iteration=3),col=3)
#' }
#' @export
vcgIsosurface <- function(vol,threshold,spacing=NULL, origin=NULL) {
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    
    mvol <- max(vol)
    minvol <- min(vol)
    if (threshold == mvol)
        threshold <- threshold-1e-5
    else if (threshold > mvol || threshold < minvol)
        stop("threshold is outside volume values")
    storage.mode(vol) <- "integer"
    volmesh <- .Call("RMarchC",vol,threshold)
    volmesh$vb <- rbind(volmesh$vb,1)
    volmesh$it <- volmesh$it
    
    class(volmesh) <- "mesh3d"
    if (!is.null(spacing))
        volmesh$vb[1:3,] <- volmesh$vb[1:3,]*spacing
    
    if (!is.null(origin))
            volmesh$vb[1:3,] <- volmesh$vb[1:3,]+origin
    
    return(volmesh)
}

