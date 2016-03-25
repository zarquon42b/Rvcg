#' Ball pivoting surface reconstruction
#'
#' Ball pivoting surface reconstruction
#' @param x k x 3 matrix or object of class mesh3d
#' @param radius The radius of the ball pivoting (rolling) over the set of points. Gaps that are larger than the ball radius will not be filled; similarly the small pits that are smaller than the ball radius will be filled. 0 = autoguess.
#' @param clustering Clustering radius (fraction of ball radius). To avoid the creation of too small triangles, if a vertex is found too close to a previous one, it is clustered/merged with it.
#' @param angle Angle threshold (radians). If we encounter a crease angle that is too large we should stop the ball rolling.
#' @param deleteFaces in case x is a mesh and \code{deleteFaces=TRUE}, existing faces will be deleted beforehand.
#' @return triangular face of class mesh3d
#' @examples
#' require(Morpho)
#' data(nose)
#' nosereko <- vcgBallPivoting(shortnose.lm)
#' @export
vcgBallPivoting <- function(x, radius=0, clustering=0.2, angle=pi/2,deleteFaces=FALSE) {
    if (is.matrix(x)) {
        x <- list(vb=t(x))
    }  else if (!inherits(x,"mesh3d")) {
        stop("x must be matrix or mesh3d")
    }
    out <- .Call("Rballpivoting",x,radius,clustering,angle,deleteFaces)
    return(out)
}
