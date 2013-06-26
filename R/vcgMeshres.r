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
#' @keywords ~kwd1 ~kwd2
#' @export vcgMeshres
vcgMeshres <- function(mesh)
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it-1
    tmp <- .Call("Rmeshres",vb,it)
    return(tmp)
  }
