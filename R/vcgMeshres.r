#' calculate average edge length of a triangular mesh
#' 
#' calculate average edge length of a triangular mesh, by iterating over all
#' faces.
#' 
#' 
#' @param mesh triangular mesh stored as object of class "mesh3d"
#' @return Item res average edge length (a.k.a. mesh resolution)
#' @return Item edgelength vector containing lengths for each edge
#' @author Stefan Schlager
#' @examples
#' data(humface)
#' mres <- vcgMeshres(humface)
#' #histogram of edgelength distribution
#' hist(mres$edgelength)
#' #visualise average edgelength
#' points( mres$res, 1000, pch=20, col=2, cex=2)
#' @keywords ~kwd1 ~kwd2
#' 
#' @export vcgMeshres
vcgMeshres <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    tmp <- .Call("Rmeshres",vb,it)
    return(tmp)
  }
