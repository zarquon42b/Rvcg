#' Subsample mesh surface
#'
#' subsample surface of a triangular mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param SampleNum integer Number of sampled points
#' @param type seclect sampling type
#' @details not ready yet
#' @return sampled points
#' @examples
#' require(rgl)
#' data(humface)
#' ss <- vcgSample(humface,SampleNum = 200)
#' points3d(ss)
#' @export vcgSample
vcgSample <- function(mesh, SampleNum=10,type=1)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        
        tmp <- .Call("Rsample", vb, it, SampleNum, type)
        tmp <- t(tmp)
        #tmp$normals <- rbind(tmp$normals,1)
        class(tmp) <- "mesh3d"
        return(tmp)
        
    }
