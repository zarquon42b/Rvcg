#' check if a mesh is intersected by a set of rays
#'
#' check if a mesh is intersected by a set of rays
#' @param x a triangular mesh of class 'mesh3d' or a list containing vertices and vertex normals (fitting the naming of 'mesh3d'.
#' @param mesh triangular mesh to be intersected.
#' @return list with following items:
#' \item{vb }{3 x n matrix containing intersection points}
#' \item{normals }{4 x n matrix containing homogenous coordinates of normals at intersection points}
#' \item{quality }{integer vector: 1 indicates that a ray has intersected 'mesh' , while 0 means not}
#' \item{distance }{numeric vector: distances to intersection}
#' 
#' @export vcgIntersect
vcgIntersect <- function(x,mesh)
{
  if (!inherits(mesh,"mesh3d") || !inherits(x,"mesh3d"))
            stop("arguments 'x' and 'mesh' needs to be object of class 'mesh3d'")
  vb <- mesh$vb[1:3,]
  it <- mesh$it - 1
  dimit <- dim(it)[2]
  dimvb <- dim(vb)[2]
  storage.mode(it) <- "integer"
  if (is.null(x$normals))
    adnormals(x)
  clost <- x$vb[1:3,]
  normals <- x$normals[1:3,]
  clostDim <- ncol(clost)
  dis <- rep(0,clostDim)
  hit <- dis
  storage.mode(hit) <- "integer"
  tmp <- .C("Rintersect",vb,ncol(vb),it,ncol(it),clost,clostDim,normals,dis,hit)
  x$vb[1:3,] <- tmp[[5]]
  x$normals <- rbind(tmp[[7]],1)
  x$quality <- tmp[[9]]
  x$distance <- tmp[[8]]
  return(x)
}
