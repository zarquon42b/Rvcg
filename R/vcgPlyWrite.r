#' Export meshes to PLY-files
#'
#' Export meshes to PLY-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d'
#' @param filename character: filename (file extension '.ply' will be added automatically.
#' @param binary logical: write binary file
#' @param addNormals logical: compute per-vertex normals and add to file
#' @param writeCol logical: export existing per-vertex color stored in mesh$material$color
#' @examples
#' data(humface)
#' vcgPlyWrite(humface,filename = "humface")
#' @export vcgPlyWrite
vcgPlyWrite <- function(mesh, filename=dataname, binary = TRUE, addNormals = FALSE, writeCol=TRUE)
{
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    it <- (mesh$it-1)
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    storage.mode(it) <- "integer"
     if ( FALSE %in% is.integer(c(it)) || FALSE %in% is.numeric(c(vb)) || !is.character(filename) )
         stop("Please provide sensible arguments!")
    filename <- paste(filename,".ply",sep="")
    if (!is.null(mesh$material$color) && writeCol==TRUE) {
        ## setup color export
        hasCol <- TRUE
        vn <- ncol(vb)
        col = rep("#FFFFFF", vn)
        tmp1 <- data.frame(it = as.vector(mesh$it))
        tmp1$rgb <- as.vector(mesh$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
        colvec <- matrix(col2rgb(col), 3, vn, byrow = F)
        storage.mode(colvec) <- "integer"
    }
    binary <- as.logical(binary)
    addNormals <- as.logical(addNormals)
    
    tmp <- .Call("RPlyWrite", vb, it , binary, addNormals, filename, colvec, hasCol)
}
