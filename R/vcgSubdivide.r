#' subdivide the triangles of a mesh
#'
#' subdivide the triangles of a mesh
#' @param x triangular mesh of class "mesh3d"
#' @param threshold minimum edge length to subdivide
#' @param type character: algorithm used. Options are Butterfly and Loop (see notes)
#' @param looptype character: method for type = loop options are "loop","regularity","continuity" (see notes)
#' @param iterations integer: number of iterations
#' @param silent logical: suppress output.
#' @return returns subdivided mesh
#' @note
#' The different algorithms are (from meshlab description):
#' \itemize{
#' \item{\bold{Butterfly Subdivision:} Apply Butterfly Subdivision Surface algorithm. It is an interpolated method, defined on arbitrary triangular meshes. The scheme is known to be C1 but not C2 on regular meshes}
#' \item{\bold{Loop Subdivision:} Apply Loop's Subdivision Surface algorithm. It is an approximant subdivision method and it works for every triangle and has rules for extraordinary vertices. Options are "loop" a simple subdivision, "regularity" to enhance the meshe's regularity and "continuity" to enhance the mesh's continuity.}
#' }
#' @examples
#' data(humface)
#' subdivide <- vcgSubdivide(humface,type="Loop",looptype="regularity")
#' 
#' @export
vcgSubdivide <- function(x, threshold=NULL, type=c("Butterfly","Loop"),looptype=c("loop","regularity","continuity"),iterations=3, silent= FALSE) {
    typeargs <- c("butterfly","loop")
    if (is.null(threshold))
        threshold <- -1
    type <- match.arg(tolower(type[1]),typeargs)
    type <- match(type,typeargs)-1
    loopargs <- c("loop","regularity","continuity")
    looptype <- match.arg(tolower(looptype[1]),loopargs)
    looptype <- match(looptype,loopargs)-1
    if (!inherits(x,"mesh3d"))
        stop("mesh must be of class mesh3d")
    out <- .Call("Rsubdivision",x,iterations,threshold,type,looptype,silent)
    return(out)
}
