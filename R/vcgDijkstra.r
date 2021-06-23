
#' @title Compute pseudo-geodesic distances on a triangular mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param vertpointer integer: references indices of vertices on the mesh
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
#' @note Make sure to have a clean manifold mesh. Note that this computes the length of the pseudo-geodesic path (following the edges) between the two vertices.
#' @export
vcgDijkstra <- function(x, vertpointer) {
    vertpointer <- as.integer(vertpointer-1)
    vb <- x$vb
    it <- x$it-1
    out <- .Call("Rdijkstra",vb,it,vertpointer)
    return(out)
}


#' @title Compute pseudo-geodesic distance between two points on a mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param pt1 3D coordinate on mesh or index of vertex
#' @param pt2 3D coordinate on mesh or index of vertex
#' @return returns the geodesic distance between \code{pt1} and \code{pt2}.
#' @note Make sure to have a clean manifold mesh. Note that this computes the length of the pseudo-geodesic path (following the edges) between the two vertices closest to these points.
#' @examples
#' data(humface)
#' pt1 <- humface.lm[1,]
#' pt2 <- humface.lm[5,]
#' vcgGeodist(humface,pt1,pt2)
#' @export
vcgGeodist <- function(x,pt1,pt2) {
    if (length(pt1) == 1)
        pt1 <- x$vb[1:3,pt1]
    if (length(pt2) == 1)
        pt2 <- x$vb[1:3,pt2]
    mypts <- rbind(pt1,pt2)
    clost <- vcgKDtree(x,mypts,k=1)
    geo <- vcgDijkstra(x,vertpointer = clost$index[1,1])[clost$index[2,1]]
    return(geo)
}


#' @title Compute geodesic path between vertices on a mesh
#' @param x triangular mesh of class \code{mesh3d}
#' @param source scalar positive integer, the source vertex index
#' @param targets positive integer vector, the target vertex indices
#' @return list of integer vectors, representing the paths.
#' @examples
#' data(humface)
#' vcgGeodesicPath(humface,1,c(5,10))
#' @export
vcgGeodesicPath <- function(x, source, targets, maxdist=1e6) {
  vertpointer_source <- as.integer(source - 1L)
  vertpointer_targets <- as.integer(targets - 1L)
  vb <- x$vb
  it <- x$it - 1L
  out <- .Call("RGeodesicPath",vb,it,vertpointer_source,vertpointer_targets, maxdist)
  return(out)
}


