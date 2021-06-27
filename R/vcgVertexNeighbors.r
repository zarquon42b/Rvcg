
#' @title Compute mesh adjacency list representation or the vertex neighborhoods of specific mesh vertices.
#'
#' @description Compute the \code{k}-ring vertex neighborhood for all query vertex indices \code{vi}. If only a mesh is passed (parameter \code{x}) and the other parameters are left at their default values, this compute the adjacency list representation of the mesh.
#'
#' @param x tmesh3d instance from the \code{rgl} package
#'
#' @param vi optional, vector of positive vertex indices for which to compute the neighborhoods. All vertices are used if left at the default value \code{NULL}.
#'
#' @param numstep positive integer, the number of times to extend the neighborhood from the source vertices (the \code{k} for computing the \code{k}-ring neighborhood). Setting this to high values significantly increases the computational cost.
#'
#' @param include_self logical, whether the returned neighborhood for a vertex \code{i} should include \code{i} itself.
#'
#' @return list of positive integer vectors, the neighborhoods.
#'
#' @examples
#' data(humface)
#' adjacency_list <- vcgVertexNeighbors(humface)
#' v500_5ring = vcgVertexNeighbors(humface, vi=c(500), numstep = 5)
#'
#' @export
vcgVertexNeighbors <- function(x, vi=NULL, numstep=1L, include_self=FALSE) {
  if(is.null(vi)) {
    vi = seq(ncol(x$vb));
  }
  vi <- as.integer(vi - 1L)
  vb <- x$vb
  it <- x$it - 1L
  out <- .Call("RVVadj",vb,it,vi,numstep,as.integer(include_self));
  return(out);
}


#' \dontrun{
#'   if(requireNamespace("fsbrain", quitely=TRUE)) {
#'   sjd = fsaverage.path(TRUE);
#'   surface = subject.surface(sjd, 'fsaverage', surface = "white", hemi = "lh");
#'   neigh = vcgVertexNeighbors(fs.surface.to.tmesh3d(surface));
#'   fsbrain::highlight.vertices.on.subject(sjd, 'fsaverage', verts_lh=neigh[[100]]);
#'   # assert max(unlist(lapply(neigh, max))) == 1643842 }
