#' calculate the Hausdorff beteween 2 meshes
#'
#' calculate the Hausdorff beteween two triangular meshes.
#'
#' @param ref triangular mesh (object of class 'mesh3d')
#' @param target triangular mesh (object of class 'mesh3d')
#'
#' @return
#' \item{maxdist}{maximal Hausdorff distance}
#' \item{meandist}{mean Hausdorff distance}
#' \item{RMSdist}{RMS of the Hausdorff distances}
#' \item{nvbsamples}{number of vertices sampled}
#' \item{nsamples}{number of samples}
#'
#' @examples
#'
#' data(humface)
#' curv <- vcgCurve(humface)
#' ##visualise per vertex mean curvature
#' \dontrun{
#' require(Morpho)
#' meshDist(humface,distvec=curv$meanvb,from=-0.2,to=0.2,tol=0.01)
#' }
#' @export vcgHausdorff
vcgHausdorff <- function(ref, target) {
   if (!inherits(ref,"mesh3d"))
       stop("argument 'ref' needs to be object of class 'mesh3d'")
   if (!inherits(target,"mesh3d"))
       stop("argument 'target' needs to be object of class 'mesh3d'")

   tmp <- .Call("Rhausdorff",ref$vb[1:3,],ref$it-1L,target$vb[1:3,],target$it-1L)
   return(tmp)
}
