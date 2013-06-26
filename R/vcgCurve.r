#' @export vcgCurve
vcgCurve <- function(mesh)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        tmp <- .Call("Rcurvature",vb,it)
        return(tmp)
  }
