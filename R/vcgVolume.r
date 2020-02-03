#' Compute volume for manifold meshes
#'
#' Compute volume for manifold meshes
#' @param x triangular mesh of class mesh3d
#' @return returns volume
#' @examples
#' mysphere <- vcgSphere()
#' vcgVolume(mysphere)
#' \dontrun{
#' ## here is an example where the mesh has some non-manifold vertices
#' 
#' mysphere <- vcgSphere(normals=FALSE)
#' ## add a degenerate face
#' mysphere$it <- cbind(mysphere$it,c(1,2,1))
#' try(vcgVolume(mysphere))
#'
#' ## fix the error using vcgClean():
#' vcgVolume(vcgClean(mysphere,sel=0:6,iterate=TRUE))
#' }
#' 
#' 
#' @note
#' Please note, that this function only works reliably on watertight, coherently oriented meshes that constitute a manifold. In case your mesh has some issues regarding non-manifoldness or there are isolated pieces flying around, you can use vcgIsolated and vcgClean to remove those.
#' @export
vcgVolume <- function(x) {
    out <-  .Call("Rmeshvol",x)
    return(out)
}

