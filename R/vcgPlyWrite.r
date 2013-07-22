#' Export meshes to PLY-files
#'
#' Export meshes to PLY-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d'
#'  @param filename character: filename (file extension '.ply' will be added automatically.
#' @param binary logical: write binary file
#' @param addNormals logical: compute per-vertex normals and add to file
#' @examples
#' data(humface)
#' vcgPlyWrite(humface,filename = "humface")
#' @export vcgPlyWrite
vcgPlyWrite <- function(mesh, filename=dataname, binary = TRUE, addNormals = FALSE)
{
    vb <- mesh$vb[1:3,]
    it <- (mesh$it-1)
    dataname <- deparse(substitute(mesh))
    filename <- as.character(filename)
    storage.mode(it) <- "integer"
     if ( FALSE %in% is.integer(c(it)) || FALSE %in% is.numeric(c(vb)) || !is.character(filename) )
         stop("Please provide sensible arguments!")
    filename <- paste(filename,".ply",sep="")
    tmp <- .Call("RPlyWrite", vb, it , binary, addNormals, filename)
}
