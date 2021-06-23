#' Isotropically remesh a triangular surface mesh
#'
#' Isotropically remesh a triangular surface mesh
#' @param x mesh of class \code{mesh3d}
#' @param TargetLen numeric: edge length of the target surface
#' @param FeatureAngleDeg define Crease angle (in degree).
#' @param MaxSurfDist Max. surface distance
#' @param iterations ToDo
#' @param Adaptive enable adaptive remeshing
#' @param split enable refine step
#' @param collapse enable collapse step
#' @param swap enable dge swap
#' @param smooth enable smoothing
#' @param project enable reprojection step
#' @param surfDistCheck check distance to surface
#' @return returns the remeshed surface mesh
#' @examples \dontrun{
#' data(humface)
#' resampledMesh <- vcgIsotropicRemeshing(humface,TargetLen=2.5)
#' }
#' @export 
vcgIsotropicRemeshing <- function(x,TargetLen=1, FeatureAngleDeg=10, MaxSurfDist=1, iterations=3, Adaptive=FALSE, split=TRUE, collapse=TRUE, swap=TRUE, smooth=TRUE,project=TRUE, surfDistCheck=TRUE) {
   
    vb <- x$vb
    it <- x$it-1
    out <- .Call("RisotropicResampling",vb,it,TargetLen, FeatureAngleDeg, MaxSurfDist, iterations, Adaptive, split, collapse, swap, smooth,project, surfDistCheck)
   
    class(out) <- "mesh3d"
    return(out)
}
