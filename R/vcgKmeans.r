#' fast Kmean clustering for 1D, 2D and 3D data
#' @param x matrix containing coordinates or mesh3d
#' @param k number of clusters
#' @param iter.max maximum number of iterations
#' @param getClosest logical: if TRUE the indices of the points closest to the k-centers are sought.
#' @param threads integer: number of threads to use
#' @return
#' returns a list containing
#' \item{centers}{cluster center}
#' \item{class}{vector with cluster association for each coordinate}
#' If \code{getClosest=TRUE}
#'  \item{selected}{vector with indices of points closest to the centers}
#' 
#' @examples
#' require(Rvcg);require(rgl)
#' data(humface)
#' set.seed(42)
#' clust <- vcgKmeans(humface,k=1000,threads=1)
#' @seealso \code{\link{vcgSample}}
#' @export
vcgKmeans <- function(x,k=10,iter.max=10,getClosest=FALSE,threads=0) {
    if (is.vector(x))
        x <- as.matrix(x)
    if (is.matrix(x)) {
        origdim <- ncol(x)
        if (origdim == 2)
            x <- cbind(x,0)
        else if (origdim == 1)
            x <- cbind(x,0,0)
        if (ncol(x) == 3) {
            x <- list(vb=t(x))
            class(x) <- "mesh3d"
        } else 
            stop("if x is a matrix, only 2 or 3 columns are allowed")
    } else {
        if(!inherits(x,"mesh3d"))
            stop("only meshes or matrices allowed")
        else
            origdim <- 3
    }
    if (k > ncol(x$vb))
        stop("number of centers exceeds sample size")
    set.seed(rnorm(1))
    out <- .Call("Rkmeans",x,k,iter.max,threads)
    
    if (getClosest)
        out$selected <- sort(unique(vcgKDtree(x,out$centers,k=1,threads = threads)$index))
    out$centers <- out$centers[,1:origdim]
    return(out)
    
}
