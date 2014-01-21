#' update vertex normals of a triangular mesh
#'
#' update vertex normals of a triangular mesh
#' @param mesh triangular mesh of class 'mesh3d'
#' @param type select the method to compute per-vertex normals: 0=area weighted average of surrounding face normals; 1 = angle weighted vertex normals.
#'
#' @return mesh with updated/created normals
#'
#' @examples
#' data(humface)
#' humface$normals <- NULL # remove normals
#' humface <- vcgUpdateNormals(humface)
#' @export

vcgUpdateNormals <- function(mesh,type = 0)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it-1
        if (! type %in% c(0,1))
            stop("please set valid type")
        normals <- .Call("RupdateNormals", vb, it, type)
        mesh$normals <- rbind(normals,1)
        return(mesh)
    }
    
