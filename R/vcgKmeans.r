#' fast Kmean clustering for 2D and 3D data
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
#'  \item{clost_center}{vector with indices of points closest to the centers}
#' 
#' @examples
#' require(Rvcg);require(rgl)
#' data(humface)
#' set.seed(42)
#' clust <- vcgKmeans(humface,k=1000,threads=2)
#' @seealso \code{\link{vcgSample}}
#' @importFrom parallel detectCores
#' @export
vcgKmeans <- function(x,k=10,iter.max=10,getClosest=FALSE,threads=parallel::detectCores()) {
    if (is.matrix(x)) {
        if (ncol(x) == 2)
            x <- cbind(x,0)
        if (ncol(x) == 3) {
            x <- list(vb=t(x))
            class(x) <- "mesh3d"
        } else 
            stop("if x is a matrix, only 2 or 3 columns are allowed")
    } else {
        if(!inherits(x,"mesh3d"))
            stop("only meshes or matrices allowed")
    }
    if (k > ncol(x$vb))
        stop("number of centers exceeds sample size")
    set.seed(rnorm(1))
    out <- .Call("Rkmeans",x,k,iter.max,threads)
    if (getClosest)
        out$clost_center <- sort(unique(vcgKDtree(x,out$centers,k=1,threads = threads)$index))
    return(out)
    
}
