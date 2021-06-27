
#' @title Compute mesh adjacency list representation or the vertex neighborhoods of specific mesh vertices.
#'
#' @param x tmesh3d instance
#'
#' @param vi optional, vector of positive vertex indices for which to compute neighborhood. All vertices are used if left at the default value \code{NULL}.
#'
#' @param numstep positive integer, the number of times to extend the neighborhood (the \code{k} for computing the k-ring neighborhood).
#'
#' @param include_self logical, whether the returned neighborhood for a vertex \code{i} should include \code{i} itself.
#'
#' @return list of integer vectors, the neighborhoods.
#'
#' @examples
#' data(humface)
#' adjacency_list <- vcgVertexNeighbors(humface)
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
