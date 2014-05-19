#' Create Isosurface from 3D-array
#'
#' Create Isosurface from 3D-array using Marching Cubes algorithm
#'
#' @param vol an integer valued 3D-array
#' @param att a list containing the additional information: a$spacing needs to be a vector of length 3 indicating the voxel dimensons.
#' @param lower numeric:lower threshold
#' @param upper numeric: upper threshold
#' @return returns a triangular mesh of class "mesh3d"
#' @examples
#' #this is the example from the package "misc3d"
#' x <- seq(-2,2,len=50)
#' g <- expand.grid(x = x, y = x, z = x)
#' v <- array(g$x^4 + g$y^4 + g$z^4, rep(length(x),3))
#' storage.mode(v) <- "integer"
#' mesh <- Rvcg:::vcgMarchingCube(v,lower=1)
#' \dontrun{
#' require(rgl)
#' wire3d(mesh)
#' }
#' @export
vcgIsosurface <- function(vol,att=NULL,lower=min(vol),upper=max(vol)) {
### vol = 3D array with numeric grey values
### att=attributes about original voxel data where att$x=spacing of x voxels, att$y, etc, and att$origin the origin of the original data.
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    #storage.mode(tol) <- "numeric"
    o <- .Call("RMarchC",vol,lower,upper)
    o$vb <- rbind(o$vb,1)
    o$it <- o$it
    storage.mode(vol) <- "integer"
    class(o) <- "mesh3d"
    if (!is.null(att)) {
        if (length(att$spacing == 3)) {
            o$vb[1,] <- o$vb[1,]*att$spacing[1]
            o$vb[2,] <- o$vb[2,]*att$spacing[2]
            o$vb[3,] <- o$vb[3,]*att$spacing[3]
        }
        if (!is.null(att$origin))
            o$vb[1:3,] <- o$vb[1:3,]+att$origin
    }
    return(o)
}

