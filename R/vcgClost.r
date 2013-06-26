#' Project coordinates onto a target triangular surface mesh.
#' 
#' For a set of 3D-coordinates the closest matches on a target surface are
#' determined and normals at as well as distances to that point are calculated.
#' 
#' 
#' @param x k x 3 matrix containing 3D-coordinates or object of class "mesh3d".
#' @param mesh triangular surface mesh stored as object of class "mesh3d".
#' @param sign logical: if TRUE, signed distances are returned.
#' @param barycentric logical: if TRUE, barycentric coordinates of the hit
#' points are returned.
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{ 4xn matrix containing n vertices as homolougous coordinates.}
#' \item{normals }{4xn matrix containing vertex normals.}
#' \item{quality }{vector: containing distances to target.}
#' \item{it }{4xm matrix containing vertex indices forming triangular
#' faces.Only available, when x is a mesh.}
#' \item{barycoords }{3xm Matrix containing barycentric coordinates of
#' closest points; only available if barycentric=TRUE.}
#' @note If large part of the reference mesh are far away from the target
#' surface, calculation can become very slow. For this case, the function
#' "closemeshKD" from the package "Morpho" is suggested.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling.
#' @keywords ~kwd1 ~kwd2
#' 
#' 
#' 
#' 
#' @export vcgClost
vcgClost <- function(x,mesh,sign=TRUE,barycentric=FALSE)
  
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it - 1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    storage.mode(it) <- "integer"

    if (is.matrix(x))
      {
        clost <- t(x)
        x <- list()
        x$vb <- clost
      }
    else
      {
        clost <- x$vb[1:3,]
      }
    
    border <- rep(0,ncol(x$vb))
    storage.mode(border) <- "integer"
    clostDim <- ncol(clost)
    faceptr <- dis <- rep(0,clostDim)
    storage.mode(clost) <- "double"
    storage.mode(faceptr) <- "integer"
    barycoord <- normals <- clost*0
    sign <- as.integer(sign)
    barycentric <- as.integer(barycentric)
    tmp <- .C("Rclost",vb,ncol(vb),it,ncol(it),clost,clostDim,clost,dis,sign,border,barycentric,barycoord,faceptr=faceptr)
    x$vb[1:3,] <- tmp[[5]]
    x$normals <- rbind(tmp[[7]],1)
    chcknorm <- which(is.nan(x$normals))
    if (length(chcknorm) > 0)
      x$normals[chcknorm] <- 0
                      
    x$quality <- tmp[[8]]
x$border <- tmp[[10]]
    if(barycentric==1)
        x$barycoords <- tmp[[12]]
    x$faceptr=tmp$faceptr+1
    invisible(x)
  }
    
