vcgMarchingCube <- function(vol,att=NULL) {
### vol = 3D array with numeric grey values
### att=attributes about original voxel data where att$x=spacing of x voxels, att$y, etc, and att$origin the origin of the original data.
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    o <- .Call("RMarchC",vol)
    o$vb <- rbind(o$vb,1)
    o$it <- o$it
    class(o) <- "mesh3d"
    if (!is.null(att)) {
        o$vb[3,] <- o$vb[3,]*diff(att$z)[1]
        o$vb[3,] <- o$vb[3,]*diff(att$z)[1]
        o$vb[3,] <- o$vb[3,]*diff(att$z)[1]
        o$vb[1:3,] <- o$vb[1:3,]+att$origin
    }
    return(o)
}

