#' Compute volume for manifold meshes
#'
#' Compute volume for manifold meshes
#' @param x triangular mesh of class mesh3d
#' @return returns volume
#' @examples
#' mysphere <- vcgSphere()
#' vcgVolume(mysphere)
#'
#' @note
#' Please note, that this function only works reliably on watertight, coherently oriented meshes that constitute a manifold. 
#' @export
vcgVolume <- function(x) {
    out <-  .Call("Rmeshvol",x)
    return(out)
}

