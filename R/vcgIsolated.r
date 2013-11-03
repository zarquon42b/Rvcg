#' Remove isolated pieces from a surface mesh.
#' 
#' Remove isolated pieces from a surface mesh, selected by face number or
#' diameter. Also the option only to keep the largest piec can be selected
#' 
#' 
#' @param mesh triangular mesh of class "mesh3d".
#' @param facenum integer: all connected pieces with less components are
#' removed. If not specified or 0 and diameter is NULL, then only the component
#' with the most faces is kept. 
#' @param diameter numeric: all connected pieces smaller diameter are removed
#' removed. 0 removes all component but the largest ones. This option overrides
#' 'facenum'.
#' @return returns the reduced mesh.
#' @author Stefan Schlager
#' @seealso \code{\link{vcgPlyRead}}
#' @keywords ~kwd1 ~kwd2
#' @examples
#' 
#' require(rgl)
#' data(humface)
#' cleanface <- vcgIsolated(humface)
#' 
#' 
#' @export vcgIsolated
vcgIsolated <- function(mesh,facenum=NULL,diameter=NULL)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        if (is.null(facenum))
            facenum <- 0
        if (!is.null(diameter))
            facenum <- -1

        if (is.null(diameter))
            diameter <- 0
        vb <- mesh$vb[1:3,]
        it <- mesh$it-1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        tmp <- .Call("Risolated", vb, it, diameter, facenum)
        outmesh <- list()
        class(outmesh) <- "mesh3d"
        outmesh$vb <- rbind(tmp$vb,1)
        outmesh$it <- tmp$it
        outmesh$normals <- rbind(tmp$normals, 1)
        invisible(outmesh)
    }
