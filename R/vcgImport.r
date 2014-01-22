#' Import common mesh file formats.
#' 
#' Import common mesh file formats and store the results in an object of
#' class "mesh3d" - momentarily only triangular meshes are supported.
#' 
#' 
#' @param file character: file to be read.
#' @param updateNormals logical: if TRUE and the imported file contais faces,
#' vertex normals will be (re)calculated. Otherwise, normals will be a matrix containing zeros.
#' @param readcolor if TRUE, vertex colors will be read if available, otherwise all vertices will be colored white.
#' @param clean if TRUE, duplicated and unreferenced vertices are removed (be careful when importing point clouds).
#' @return Object of class "mesh3d"
#' 
#' with:
#' \item{vb }{4 x n matrix containing n vertices as homolougous coordinates}
#' \item{it }{3 x m matrix containing vertex indices forming triangular faces}
#' \item{normals }{4 x n matrix containing vertex normals (homologous coordinates)}
#' @author Stefan Schlager
#' @seealso \code{\link{vcgSmooth}}
#' @examples
#' data(humface)
#' vcgPlyWrite(humface)
#' readit <- vcgImport("humface.ply")
#' @keywords ~kwd1 ~kwd2
#' @export 

vcgImport <- function(file, updateNormals = TRUE, readcolor=FALSE, clean = TRUE) {
    ncfile <- nchar(file)
    ext <- substr(file,ncfile-2,ncfile)
    file <- path.expand(file)
    x <- file
    if (! file.exists(x))
        stop("no such file")
    updateNormals <- as.logical(updateNormals)
    readcolor <- as.logical(readcolor)
    clean <- as.logical(clean)


    tmp <- .Call("RallRead", file, updateNormals, readcolor, clean)
    if (!is.list(tmp))
        stop("mesh is not readable")
    out <- list()
    class(out) <- "mesh3d"
    out$vb <- rbind(matrix(tmp$vb,3,length(tmp$vb)/3),1)
    out$it <- matrix(tmp$it,3,(length(tmp$it)/3))+1
    out$normals <- rbind(matrix(tmp$normals,3,length(tmp$normals)/3),1)
    if (readcolor)
        {
          colvec <- out$colors
          out$material <- list()
          colvec <- rgb(colvec[1,],colvec[2,],colvec[3,],maxColorValue=255)
          colfun <- function(x)
            {
              x <- colvec[x]
              return(x)
            }
          out$material$color <- matrix(colfun(out$it),dim(out$it))
        }
    return(out)
}
