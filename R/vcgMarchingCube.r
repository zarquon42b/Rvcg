vcgMarchingCube <- function(vol,att=NULL,tol=1) {
### vol = 3D array with numeric grey values
### att=attributes about original voxel data where att$x=spacing of x voxels, att$y, etc, and att$origin the origin of the original data.
    if (length(dim(vol)) != 3)
        stop("3D array needed")
    storage.mode(tol) <- "numeric"
    o <- .Call("RMarchC",vol,tol=tol)
    o$vb <- rbind(o$vb,1)
    o$it <- o$it
    storage.mode(vol) <- "integer"
    class(o) <- "mesh3d"
    if (!is.null(att)) {
        o$vb[1,] <- o$vb[1,]*diff(att$x)[1]
        o$vb[2,] <- o$vb[2,]*diff(att$y)[1]
        o$vb[3,] <- o$vb[3,]*diff(att$z)[1]
        if (!is.null(att$origin))
            o$vb[1:3,] <- o$vb[1:3,]+att$origin
    }
    return(o)
}

