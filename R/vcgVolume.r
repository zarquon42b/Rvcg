#' Compute volume for manifold meshes
#'
#' Compute volume for manifold meshes
#' @param x triangular mesh of class mesh3d
#' @return returns volume
#' @examples
#' mysphere <- vcgSphere()
#' vcgVolume(mysphere)
#' 
#' @export
vcgVolume <- function(x) {
    out <-  .Call("Rmeshvol",x)
    return(out)
}

