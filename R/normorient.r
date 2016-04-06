cSizeMesh <- function(mesh) {
    x <- t(mesh$vb[1:3,])
    X <- scale(x, scale = FALSE)
    y <- sqrt(sum(as.vector(X)^2))
    return(y)
}
meshOff <- function(x,offset) {
    x <- vcgUpdateNormals(x)
    x$vb[1:3,] <- x$vb[1:3,]+offset*x$normals[1:3,]
    return(x)
}
#' check the orientation of a mesh
#'
#' check the orientation of a mesh assuming that expansion along normals increases centroid size
#'
#' @param x mesh of class mesh3d
#' @param offset numeric: amount to offset the mesh along the vertex normals. If NULL a reasonable value will be estimated.
#' @return returns TRUE if mesh is oriented correctly and FALSE otherwise
#' @details assuming that a correctly (i.e outward) oriented mesh increases its centroid size when 'growing' outwards, this function tests whether this is the case.
#' @examples
#' data(dummyhead)
#' ## now we invert faces inwards
#' checkFaceOrientation(dummyhead.mesh)
#' dummyinward <- Morpho::invertFaces(dummyhead.mesh)
#' checkFaceOrientation(dummyinward)
#' @export
checkFaceOrientation <- function(x,offset=NULL) {
    if (is.null(offset))
        offset <- pcax(x)/15
    out <- TRUE
    xoff <- meshOff(x,offset)
    cx <- cSizeMesh(x)
    cxoff <- cSizeMesh(xoff)
    if (cx > cxoff)
        out <- FALSE
    return(out)
}

pcax <- function(mesh) {
    x <- t(mesh$vb[1:3,])
    pc <- 3*prcomp(x,retx=FALSE)$sdev[1]
    return(pc)
}
    
