#' get all faces belonging to each vertex in a mesh
#'
#' get all faces belonging to each vertex in a mesh
#' @param mesh triangular mesh
#'
#' @return list containing one vector per vertex with the indices of the adjacent faces 
#' @export vertAdj
#' 
vertAdj <- function(mesh) {
    it <- mesh$it-1
    storage.mode(it) <- "integer"
    vn <- ncol(mesh$vb)
    out <- .Call("RvertAdj",vn,it)
    
    return(out)
}
