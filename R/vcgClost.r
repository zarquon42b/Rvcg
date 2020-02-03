#' Project coordinates onto a target triangular surface mesh.
#' 
#' For a set of 3D-coordinates/triangular mesh, the closest matches on a
#' target surface are determined and normals at as well as distances to
#' that point are calculated.
#' 
#' 
#' @param x k x 3 matrix containing 3D-coordinates or object of class "mesh3d".
#' @param mesh triangular surface mesh stored as object of class "mesh3d".
#' @param sign logical: if TRUE, signed distances are returned.
#' @param barycentric logical: if TRUE, barycentric coordinates of the hit
#' points are returned.
#' @param smoothNormals logical: if TRUE, laplacian smoothed normals are used.
#' @param borderchk logical: request checking if the hit face is at the border of the mesh.
#' @param tol maximum distance to search. If distance is beyond that, the original point will be kept and the distance set to NaN. If tol = 0, tol is set to 2*diagonal of the bounding box of \code{mesh}.
#' @param facenormals logical: if TRUE only the facenormal of the face the closest point has hit is returned, the weighted average of the surrounding vertex normals otherwise.
#' @param ... additional parameters, currently unused.
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{4 x n matrix containing n vertices as homolougous coordinates.}
#' \item{normals }{4 x n matrix containing vertex normals.}
#' \item{quality }{numeric vector containing distances to target.}
#' \item{it }{3 x m integer matrix containing vertex indices forming triangular
#' faces.Only available, when x is a mesh.}
#' \item{border }{integer vector of length n: if borderchk = TRUE, for each clostest point the value will be 1 if the hit face is at the border of the target mesh and 0 otherwise.} 
#' \item{barycoords }{3 x m Matrix containing barycentric coordinates of
#' closest points; only available if barycentric=TRUE.}
#' \item{faceptr }{vector of face indeces on which the closest points are located}
#' @note If large part of the reference mesh are far away from the target
#' surface, calculation can become very slow. In that case, the function
#' \code{vcgClostKD} will be significantly faster.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @references Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed
#' Distance Fields From Triangle Meshes. Informatics and Mathematical
#' Modelling.
#'
#' @examples
#' data(humface)
#' clost <- vcgClost(humface.lm, humface)
#' 
#' 
#' 
#' 
#' @export vcgClost
vcgClost <- function(x,mesh,sign=TRUE,barycentric=FALSE, smoothNormals=FALSE, borderchk = FALSE,tol=0,facenormals=FALSE,...)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        mesh <- meshintegrity(mesh,facecheck=TRUE)
                
        if (is.matrix(x)) {
            clost <- t(x[,1:3,drop=FALSE])
            x <- list()
            x$vb <- clost
            class(x) <- "mesh3d"
        } else if (inherits(x,"mesh3d")) {
            clost <- x$vb[1:3,,drop=FALSE]
            if (!is.matrix(clost))
            stop("x has no vertices")
        } else
            stop("x must be a matrix or an object of class mesh3d")

        if (FALSE %in% is.logical(c(sign, barycentric, smoothNormals,borderchk)))
            stop("please provide sensible input")
        
        storage.mode(clost) <- "double"
        outmesh <- .Call("Rclost",mesh, clost,sign,borderchk,barycentric,smoothNormals,tol,facenormals)
        
        chcknorm <- which(is.nan(outmesh$normals))
        if (length(chcknorm) > 0)
            outmesh$normals[chcknorm] <- 0
        outmesh$it <- x$it
        nancheck <- which(is.nan(outmesh$quality))
        if (length(nancheck))
            warning("some points where beyond the search tolerance")
        invisible(outmesh)
    }

