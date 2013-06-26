#' Clean triangular surface meshes
#'
#' Apply several cleaning algorithms to surface meshes
#' @param mesh triangular object of class 'mesh3d
#' @param sel integer select cleaning type (see "details"
#' @details available options are 0= only duplicated vertices and faces are removed (always applied before cleaning). 1=todo.
#' @return cleaned mesh
#' @examples
#' data(humface)
#' cleanface <- vcgClean(humface, sel=1)
#' @export vcgClean
vcgClean <- function(mesh, sel = 0)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        storage.mode(sel) <- "integer"
        tmp <- .Call("Rclean", vb, it, sel)
        tmp$vb <- rbind(tmp$vb,1)
        tmp$normals <- rbind(tmp$normals,1)
        class(tmp) <- "mesh3d"
        return(tmp)
        
    }
