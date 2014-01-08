#' check if a mesh is intersected by a set of rays
#'
#' check if a mesh is intersected by a set of rays
#' @param x a triangular mesh of class 'mesh3d' or a list containing vertices and vertex normals (fitting the naming of 'mesh3d'.
#' @param mesh triangular mesh to be intersected.
#' @param tol minimum distance to target mesh
#' @param mindist search both ways (ray and -ray) and select closest point.
#' @details project a mesh along a set of given rays (stored as normals) onto a target and return the hit points as well as information if the target mesh was hit at all. If nothing is hit along the ray, the original point's value will be retrurned. If the point is already on the surface distance will be 0 but be considered as non intersecting along its normal.
#' @return list with following items:
#' \item{vb }{4 x n matrix containing intersection points}
#' \item{normals }{4 x n matrix containing homogenous coordinates of normals at intersection points}
#' \item{quality }{integer vector containing a value for each vertex of \code{x}: 1 indicates that a ray has intersected 'mesh' , while 0 means not}
#' \item{distance }{numeric vector: distances to intersection}
#' @examples
#' data(humface)
#' #get normals of landmarks
#' lms <- vcgClost(humface.lm, humface)
#' # offset landmarks along their normals for a negative amount of -5mm
#' lms$vb[1:3,] <- lms$vb[1:3,]+lms$normals[1:3,]*-5
#' intersect <- vcgRaySearch(lms, humface)
#' \dontrun{
#' require(Morpho)
#' require(rgl)
#' spheres3d(vert2points(lms),radius=0.5,col=3)
#' plotNormals(lms,long=5)
#' spheres3d(vert2points(intersect),col=2) #plot intersections
#' wire3d(humface,col="white")#'
#' }
#'
#' @export vcgRaySearch
vcgRaySearch <- function(x, mesh, tol=0, mindist=FALSE)
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
  tmp <- .Call("Rintersect",vb,it,clost,normals,tol,mindist)
  x$vb <- rbind(tmp$vb,1)
  x$normals <- rbind(tmp$normals,1)
  x$quality <- tmp$hitbool
  x$distance <- tmp$dis
  return(x)
}
