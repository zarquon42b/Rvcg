
#' create a KD-tree
#'
#' create a KD-tree
#' @param mesh matrix or triangular mesh containing coordinates
#' @param nofPointsPerCell number of points per kd-cell
#' @param maxDepth maximum tree depth
#' @return
#' returns an object of class vcgKDtree containing external pointers to the tree and the target points
#' @examples
#' data(humface)
#' mytree <- vcgCreateKDtree(humface)
#' @seealso \code{\link{vcgSearchKDtree}}
#' @export
vcgCreateKDtree <- function(mesh, nofPointsPerCell=16,maxDepth=64) {
    if (is.matrix(mesh)) {
       if (ncol(mesh) == 2)
            mesh <- cbind(mesh,0)
        if (ncol(mesh) == 3) {
            mesh <- list(vb=t(mesh))
            class(mesh) <- "mesh3d"
        } else 
            stop("if query is a matrix, only 2 or 3 columns are allowed")
       
    }
    out <- .Call("createKDtree",mesh,nofPointsPerCell,maxDepth)
    class(out) <- "vcgKDtree"
    return(out)
}

#' search an existing KD-tree
#'
#' search an existing KD-tree
#' @param kdtree object of class vcgKDtree
#' @param query atrix or triangular mesh containing coordinates
#' @param k number of k-closest neighbours to query
#' @param threads integer: number of threads to use
#' @return
#' a list with
#' \item{index}{integer matrices with indeces of closest points}
#' \item{distances}{corresponding distances}
#' @examples
#' \dontrun{
#' data(humface);data(dummyhead)
#' mytree <- vcgCreateKDtree(humface)
#' ## get indices and distances for 10 closest points.
#' closest <- vcgSearchKDtree(mytree,dummyhead.mesh,k=10,threads=1)
#' }
#' @seealso \code{\link{vcgCreateKDtree}}
#' @export
vcgSearchKDtree <- function(kdtree, query,k ,threads=0) {
    if (!inherits(kdtree,"vcgKDtree"))
        stop("no valid kdtree")
    if (is.matrix(query)) {
        if (ncol(query) == 2)
            query <- cbind(query,0)
        if (ncol(query) == 3) {
            query <- list(vb=t(query))
            class(query) <- "mesh3d"
        } else 
            stop("if query is a matrix, only 2 or 3 columns are allowed")
    } else {
        if(!inherits(query,"mesh3d"))
            stop("only meshes or matrices allowed")
    }
    out <- .Call("RsearchKDtree",kdtree$kdtree,kdtree$target,query,k,threads)
    out$index <- out$index+1
    return(out)
}

#' create a KD-tree from Barycenters for multiple closest point searches on a mesh
#'
#' create a KD-tree from Barycenters for multiple closest point searches on a mesh
#' @param mesh matrix or triangular mesh containing coordinates
#' @param nofPointsPerCell number of points per kd-cell
#' @param maxDepth maximum tree depth
#' @return
#' returns an object of class vcgKDtreeWithBarycenters containing external pointers to the tree, the barycenters and the target mesh
#' @seealso \code{\link{vcgClostOnKDtreeFromBarycenters}, \link{vcgSearchKDtree},  \link{vcgCreateKDtree}}
#' @examples
#' \dontrun{
#' data(humface);data(dummyhead)
#' barytree <- vcgCreateKDtreeFromBarycenters(humface)
#' closest <- vcgClostOnKDtreeFromBarycenters(barytree,dummyhead.mesh,k=50,threads=1)
#' }
#' @export
vcgCreateKDtreeFromBarycenters <- function(mesh, nofPointsPerCell=16,maxDepth=64) {
    mesh <- meshintegrity(mesh,facecheck=TRUE)
    barycenters <- vcgBary(mesh)
    out <- vcgCreateKDtree(barycenters,nofPointsPerCell,maxDepth)
    out$targetptr <- .Call("RmeshXPtr",mesh)
    class(out) <- "vcgKDtreeWithBarycenters"
    return(out)
}

#' search a KD-tree from Barycenters for multiple closest point searches on a mesh
#'
#' search a KD-tree from Barycenters for multiple closest point searches on a mesh
#' @param x object of class "vcgKDtreeWithBarycenters"
#' @param query matrix or triangular mesh containing coordinates
#' @param k integer: check the kdtree for the\code{k} closest faces (using faces' barycenters).
#' @param sign logical: if TRUE, signed distances are returned.
#' @param barycentric logical: if TRUE, barycentric coordinates of the hit
#' points are returned.
#' @param borderchk logical: request checking if the hit face is at the border of the mesh.
#' @param angdev maximum deviation between reference and target normals. If the none of the k closest triangles match this criterion, the closest point on the closest triangle is returned but the corresponding distance in $quality is set to 1e5.
#' @param weightnorm logical if angdev is set, this requests the normal of the closest points to be estimated by weighting the surrounding vertex normals. Otherwise, simply the hit face's normal is used (faster but slightly less accurate)
#' @param facenormals logical: if TRUE only the facenormal of the face the closest point has hit is returned, the weighted average of the surrounding vertex normals otherwise.
#' @param threads integer: threads to use in closest point search.
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{4 x n matrix containing n vertices as homolougous coordinates.}
#' \item{normals }{4 x n matrix containing vertex normals.}
#' \item{quality }{numeric vector containing distances to target.}
#' \item{it }{3 x m integer matrix containing vertex indices forming triangular
#' faces.Only available, when x is a mesh.}
#' \item{border }{integer vector of length n: if borderchk = TRUE, for each clostest point the value will be 1 if the hit face is at the border of the target mesh and 0 otherwise.} 
#' \item{barycoords }{3 x m Matrix containing barycentric coordinates of
#' closest points; only available if barycentric=TRUE.}
#' @author Stefan Schlager
#' @seealso \code{\link{vcgCreateKDtreeFromBarycenters}, \link{vcgSearchKDtree},  \link{vcgCreateKDtree}}
#' @export
vcgClostOnKDtreeFromBarycenters <- function(x,query,k=50,sign=TRUE,barycentric=FALSE, borderchk = FALSE, angdev=NULL, weightnorm=FALSE, facenormals=FALSE,threads=1) {
    if (! inherits(x,"vcgKDtreeWithBarycenters"))
        stop("provide valid object")
    if (is.matrix(query) && is.numeric(query)) {
        query <- list(vb=t(query))
        class(query) <- "mesh3d"
    } else if (!inherits(query,"mesh3d")) {
        stop("argument 'mesh' needs to be object of class 'mesh3d'")
    }
    if (is.null(angdev) || is.null(query$it))
        angdev <- 0
    out <- .Call("RsearchKDtreeForClosestPoints",x$kdtree,x$target,x$targetptr,query,k,sign,borderchk,barycentric,angdev,weightnorm,facenormals,threads)
    out$it <- query$it
    return(out)
}

