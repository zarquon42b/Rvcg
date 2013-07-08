#' Subsample mesh surface
#'
#' subsample surface of a triangular mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param SampleNum integer Number of sampled points
#' @param type seclect sampling type (1=MonteCarlo Sampling, 2=PoissonDisk Sampling)
#' @details not ready yet
#' @return sampled points
#' @examples
#' require(rgl)
#' data(humface)
#' ss <- vcgSample(humface,SampleNum = 500, type=2)
#' points3d(ss)
#' @export vcgSample
vcgSample <- function(mesh, SampleNum=10,type=1,MCsamp=20,geodes=TRUE)
    {
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
        
        return(tmp)
    }
