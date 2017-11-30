#' compute surface area of a triangular mesh
#'
#' compute surface area of a triangular mesh
#' @param mesh triangular mesh of class mesh3d
#' @param perface logical: if TRUE, a list containing the overall area, as well as the individual per-face area are reported.
#' @return surface area of mesh
#' @examples
#' data(humface)
#' vcgArea(humface)
#' @export
vcgArea <- function(mesh,perface=FALSE) {
    out <- .Call("Rarea",mesh,perface)
    return(out)
}
