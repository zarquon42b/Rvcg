#' compute surface area of a triangular mesh
#'
#' compute surface area of a triangular mesh
#' @param mesh triangular mesh of class mesh3d
#' @return surface area of mesh
#' @examples
#' data(humface)
#' vcgArea(humface)
#' @export
vcgArea <- function(mesh) {
    out <- .Call("Rarea",mesh)
    return(out)
}
