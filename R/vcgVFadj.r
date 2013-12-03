#' get all faces belonging to each vertex in a mesh
#'
#' get all faces belonging to each vertex in a mesh
#' @param mesh triangular mesh
#'
#' @return list containing one vector per vertex with the indices of the adjacent faces 
#' @export vcgVFadj
#' 
vcgVFadj <- function(mesh) {
    it <- mesh$it-1
    storage.mode(it) <- "integer"
    vb <- mesh$vb
    out <- .Call("RVFadj",vb,it)
    
    return(out)
}
