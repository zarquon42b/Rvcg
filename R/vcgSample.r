#' Subsamples points on a mesh surface
#'
#' Subsamples surface of a triangular mesh and returns a set of points located on that mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param SampleNum integer: number of sampled points (see \code{details} below)
#' @param type character: seclect sampling type ("mc"=MonteCarlo Sampling, "pd"=PoissonDisk Sampling,"km"=kmean clustering)
#' @param MCsamp integer: MonteCarlo sample iterations used in PoissonDisk sampling.
#' @param geodes logical: maximise geodesic distance between sample points (only for Poisson Disk sampling)
#' @param strict logical: if \code{type="pd"} and the amount of coordinates exceeds \code{SampleNum},  the resulting coordinates will be subsampled again by kmean clustering to reach the requested number.
#' @param iter.max integer: maximum iterations to use in k-means clustering.
#' @param threads integer number of threads to use for k-means clustering
#' @details Poisson disk subsampling will not generate the exact amount of coordinates specified in \code{SampleNum}, depending on \code{MCsamp} the result will contain more or less coordinates.
#' @return sampled points
#' @examples
#' 
#' data(humface)
#' ss <- vcgSample(humface,SampleNum = 500, type="km",threads=1)
#' \dontrun{
#' require(rgl)
#' points3d(ss)
#' }
#' @export vcgSample
vcgSample <- function(mesh, SampleNum=100,type=c("km","pd","mc"),MCsamp=20,geodes=TRUE,strict=FALSE,iter.max=100,threads=0)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        type <- match.arg(type[1],c("km","pd","mc"))
        if (type == "mc")
            type <- 1
        else if (type == "pd")
            type <- 2
        else if (type == "km")
            type <- 3
        noit <- FALSE
        if (!is.matrix(mesh$it)) {
            warning("mesh has no faces, kmean clustering on vertices is used")
            noit <- TRUE
            type <- 3
        }
        if (type %in% 1:2) {
            mesh <- meshintegrity(mesh)
            type <- as.integer(type)
            tmp <- .Call("Rsample", mesh, SampleNum, type, MCsamp, geodes)
            tmp <- t(tmp)
            if (strict && nrow(tmp) > SampleNum) {
                tmp <- vcgKmeans(tmp,k=SampleNum, iter.max=iter.max,threads=threads)$centers
                t(vcgClostKD(tmp, mesh,sign=FALSE,threads = threads)$vb[1:3,])
            }
        } else {
            tmp <-  vcgKmeans(mesh,k=SampleNum, iter.max=iter.max,threads=threads)$centers
            if (!noit)
                tmp <- t(vcgClostKD(tmp, mesh,sign=FALSE,threads=threads)$vb[1:3,])
        }
        return(tmp)
    }
