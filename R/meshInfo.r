#' get number of vertices from a mesh
#'
#' get number of vertices from a mesh
#' @param x triangular mesh
#' @return integer: number of vertices
#' @export
nverts <- function(x) {
    if (!inherits(x,"mesh3d"))
        stop("x must be of class mesh3d")
    return(ncol(x$vb))
}

#' get number of vertices from a mesh
#'
#' get number of vertices from a mesh
#' @param x triangular mesh
#' @return integer: number of triangular faces
#' @export
nfaces <- function(x) {
    if (!inherits(x,"mesh3d"))
        stop("x must be of class mesh3d")
    return(ncol(x$it))
}

#' print number of vertices and triangular faces of a mesh
#'
#' print number of vertices and triangular faces of a mesh
#' @param x triangular mesh
#' @export
meshInfo <- function(x) {
    if (!inherits(x,"mesh3d"))
        stop("x must be of class mesh3d")
    cat(paste0("mesh has ",prettyNum(nverts(x),big.mark = ",")," vertices and ",prettyNum(nfaces(x),big.mark = ",")," triangular faces\n"))
}
