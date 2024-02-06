#' Smoothes a triangular mesh
#'
#' Applies different smoothing algorithms on a triangular mesh.
#'
#' The algorithms available are Taubin smoothing, Laplacian smoothing and an
#' improved version of Laplacian smoothing ("HClaplace"). Also available are
#' Scale dependent laplacian smoothing ("fujiLaplace") and Laplacian angle
#' weighted smoothing ("angWeight")
#'
#' @param mesh triangular mesh stored as object of class "mesh3d".
#' @param type character: select smoothing algorithm. Available are "taubin",
#' "laplace", "HClaplace", "fujiLaplace", "angWeight" (and any sensible
#' abbreviations).
#' @param iteration integer: number of iterations to run.
#' @param lambda numeric: parameter for Taubin smooth (see reference below).
#' @param mu numeric:parameter for Taubin smooth (see reference below).
#' @param delta numeric: parameter for Scale dependent laplacian smoothing (see
#' reference below).and maximum allowed angle (in radians) for deviation between normals Laplacian (surface preserving).
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
#' \item{normals}{4xn matrix containing vertex normals.}
#' \item{quality }{vector: containing distances to target.}
#' \item{it }{4xm matrix containing vertex indices forming triangular
#' faces.}
#' @note The additional parameters for taubin smooth are hardcoded to the
#' default values of meshlab, as they appear to be the least distorting
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead},\link{vcgClean}}
#' @references Taubin G. 1995. Curve and surface smoothing without shrinkage.
#' In Computer Vision, 1995. Proceedings., Fifth International Conference on,
#' pages 852 - 857.
#'
#' Vollmer J., Mencl R. and Mueller H. 1999. Improved Laplacian Smoothing of
#' Noisy Surface Meshes. Computer Graphics Forum, 18(3):131 - 138.
#'
#' Schroeder, P. and Barr, A. H. (1999). Implicit fairing of irregular meshes
#' using diffusion and curvature flow: 317-324.
#' @examples
#'
#' data(humface)
#' smoothface <- vcgSmooth(humface)
#' ## view
#' \dontrun{
#' require(rgl)
#' shade3d(smoothface, col=3)
#' }
#'
#' @export vcgSmooth
vcgSmooth <- function(mesh,type=c("taubin","laplace","HClaplace","fujiLaplace","angWeight","surfPreserveLaplace"),iteration=10,lambda=0.5,mu=-0.53,delta=0.1)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        mesh <- meshintegrity(mesh)
        type <- substring(type[1],1L,1L)
        vb <- mesh$vb[1:3,,drop=FALSE]
        it <- (mesh$it-1)
        storage.mode(it) <- "integer"
        method <- 0
        if (type == "l" || type == "L") {
            method <- 1
        } else if (type == "H" || type == "h") {
            method <- 2
        } else if (type == "f" || type == "F") {
            method <- 3
        } else if (type == "a" || type == "A") {
            method <- 4
        } else if (type == "s" || type == "S") {
            method <- 5
        }
        iter <- as.integer(iteration)
        if ( FALSE %in% is.integer(c(it,iter)) || FALSE %in% is.numeric(c(vb,lambda,mu,delta)))
            stop("Please provide sensible arguments!")

        tmp <- .Call("Rsmooth",vb,it,iter,method,lambda,mu,delta)
        mesh$vb[1:3,] <- tmp$vb
        mesh$normals <- rbind(tmp$normals, 1)
        mesh$it <- tmp$it
        invisible(meshintegrity(mesh))
    }


#' Implicit Smoothes a triangular mesh
#'
#' Applies implicit smoothing algorithms on a triangular mesh.
#'
#' @param mesh triangular mesh stored as object of class "mesh3d".
#' @param lambda numeric: the amount of smoothness, useful only if
#' \code{useMassMatrix} is \code{TRUE}; default is \code{0.2}
#' @param useMassMatrix logical: whether to use mass matrix to keep the mesh
#' close to its original position (weighted per area distributed on vertices);
#' default is \code{TRUE}
#' @param fixBorder logical: whether to fix the border vertices of the mesh;
#' default is \code{FALSE}
#' @param useCotWeight logical: whether to use cotangent weight; default is
#' \code{FALSE} (using uniform 'Laplacian')
#' @param lapWeight numeric: weight when \code{useCotWeight} is \code{FALSE};
#' default is \code{1.0}
#' @param degree integer: degrees of 'Laplacian'; default is \code{1}
#' @param SmoothQ logical: whether to smooth the quality (distances to target).
#' @return returns an object of class "mesh3d" with:
#' \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
#' \item{normals}{4xn matrix containing vertex normals.}
#' \item{it }{4xm matrix containing vertex indices forming triangular
#' faces.}
#' @author Zhengjia Wang
#' @seealso \code{\link{vcgPlyRead},\link{vcgClean},\link{vcgSmooth}}
#' @examples
#'
#' data(humface)
#' smoothface <- vcgSmoothImplicit(humface)
#' ## view
#' \dontrun{
#' require(rgl)
#' shade3d(smoothface, col=3)
#' }
#'
#' @export vcgSmoothImplicit
vcgSmoothImplicit <- function(mesh,lambda=0.2,useMassMatrix=TRUE,fixBorder=FALSE,useCotWeight=FALSE,degree=1L,lapWeight=1.0,SmoothQ=FALSE)
{
    if (!inherits(mesh,"mesh3d"))
        stop("argument 'mesh' needs to be object of class 'mesh3d'")
    mesh <- meshintegrity(mesh)
    vb <- mesh$vb[1:3,,drop=FALSE]
    it <- (mesh$it-1)
    storage.mode(it) <- "integer"
    lambda <- as.double(lambda)[[1]]
    lapWeight <- as.double(lapWeight)[[1]]
    useMassMatrix <- as.logical(useMassMatrix)[[1]]
    fixBorder <- as.logical(fixBorder)[[1]]
    useCotWeight <- as.logical(useCotWeight)[[1]]
    SmoothQ <- as.logical(SmoothQ)[[1]]
    degree <- as.integer(degree)[[1]]
    tmp <- .Call("RsmoothImplicit",vb,it,lambda,useMassMatrix,fixBorder,useCotWeight,degree,lapWeight,SmoothQ)
    mesh$vb[1:3,] <- tmp$vb
    mesh$normals <- rbind(tmp$normals, 1)
    mesh$it <- tmp$it
    invisible(meshintegrity(mesh))
}

