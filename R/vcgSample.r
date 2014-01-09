#' Subsample mesh surface
#'
#' subsample surface of a triangular mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param SampleNum integer Number of sampled points
#' @param type seclect sampling type ("mc"=MonteCarlo Sampling, "pd"=PoissonDisk Sampling,"km"=kmean clustering)
#' @param MCsamp MonteCarlo sample iterations used in PoissonDisk sampling.
#' @param geodes maximise geodesic distance between sample points (only for Poisson Disk sampling)
#' @param strict if \code{type="pd"} and the amount of coordinates exceeds \code{SampleNum},  the resulting coordinates will be subsampled again by kmean clustering to reach the requested number.
#' @details Poisson disk subsampling will not generate the exact amount of coordinates specified in \code{SampleNum}, depending on \code{MCsamp} the result wil bee more or less coordinates.
#' @return sampled points
#' @examples
#' 
#' data(humface)
#' ss <- vcgSample(humface,SampleNum = 500, type="pd")
#' \dontrun{
#' require(rgl)
#' points3d(ss)
#' }
#' @export vcgSample
vcgSample <- function(mesh, SampleNum=100,type=c("km","pd","mc"),MCsamp=20,geodes=TRUE,strict=FALSE)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        type <- type[1]
        if (type == "mc")
            type <- 1
        else if (type == "pd")
            type <- 2
        else if (type == "km")
            type <- 3

        if (type %in% 1:2) {
            
            vb <- mesh$vb[1:3,]
            it <- mesh$it - 1
            dimit <- dim(it)[2]
            dimvb <- dim(vb)[2]
            storage.mode(it) <- "integer"
            type <- as.integer(type)
            SampleNum <- as.integer(SampleNum)
            MCsamp <- as.integer(MCsamp)
            if (!is.logical(geodes) || (FALSE %in% is.integer(c(it,type, MCsamp, SampleNum))) || (FALSE %in% is.numeric(vb)))
                stop("Please provide sensible arguments!")
            tmp <- .Call("Rsample", vb, it, SampleNum, type, MCsamp, geodes)
            tmp <- t(tmp)
            if (strict && nrow(tmp) > SampleNum)
                tmp <- kmeans(tmp,centers=SampleNum, iter.max=100)$centers

        } else {
            tmp <- kmeans(t(mesh$vb[1:3,]),centers=SampleNum, iter.max=100)$centers
            tmp <- t(vcgClost(tmp, mesh)$vb[1:3,])
        }
        return(tmp)
    }
