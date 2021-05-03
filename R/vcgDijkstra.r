#' Compute geodesic distances on a triangular mesh
#' 
#' Compute geodesic distances on a triangular mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param vertpointer integer: references indices of vertices on the mesh
#' @param tol numeric: threshold for max distances to consider
#' @return returns a vector of shortest distances for each of the vertices to one of the vertices referenced in \code{vertpointer}
#' @examples
#' ## Compute geodesic distance between all mesh vertices and the first vertex of a mesh
#' data(humface)
#' humface <- vcgIsolated(vcgClean(humface,sel=0:6,iterate=TRUE))
#' geo <- vcgDijkstra(humface,1)
#' if (interactive()) {
#' require(Morpho);require(rgl)
#' meshDist(humface,distvec = geo)
#' spheres3d(vert2points(humface)[1,],col=2)
#' }
#' @note Make sure to have a clean manifold mesh.
#' @export
vcgDijkstra <- function(x, vertpointer,tol=1e6) {
    vertpointer <- as.integer(vertpointer-1)
    vb <- x$vb
    it <- x$it-1
    out <- .Call("Rdijkstra",vb,it,vertpointer,tol)
    return(out)
}
