#' Resample a mesh uniformly
#'
#' Resample a mesh uniformly
#' 
#' @param x triangular mesh
#' @param voxelSize voxel size for space discretization
#' @param offset see meshlab
#' @param discretize see meshlab
#' @param multiSample see meshlab
#' @param absDist see meshlab
#' @param mergeClost logical: merge close vertices
#' @param silent logical: suppress messages
#' @return resampled mesh
#'
#' @export
vcgUniformRemesh <- function(x,voxelSize=NULL,offset=0, discretize=FALSE, multiSample=FALSE,absDist=FALSE, mergeClost=FALSE,silent=FALSE) {
    if (is.null(voxelSize))
        voxelSize <- bbox(x)$dia/50
    vb <- x$vb
    it <- x$it-1
    out <- .Call("RuniformResampling",vb,it,voxelSize,offset,discretize,multiSample, absDist,mergeClost,silent)
    out$vb <- rbind(out$vb,1)
    out$normals <- rbind(out$normals,1)
    class(out) <- "mesh3d"
    return(out)
}
bbox <- function(x) {
    bbox <- apply(t(x$vb[1:3,]), 2, range)
    bbox <- expand.grid(bbox[, 1], bbox[, 2], bbox[, 3])
    dia <- max(dist(bbox))
    return(list(bbox=bbox,diag=dia))
}
    
