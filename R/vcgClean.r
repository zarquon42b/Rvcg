#' Clean triangular surface meshes
#'
#' Apply several cleaning algorithms to surface meshes
#' @param mesh triangular mesh of class 'mesh3d'
#' @param sel integer select cleaning type (see "details")
#' @param tol numeric value determining Vertex Displacement Ratio used for splitting non-manifold vertices.
#' @details available options are: 0 = only duplicated vertices and faces are removed (always applied before cleaning). 1=remove unreferenced vertices, 2 = Remove non-manifold Faces, 3 = Remove degenerate faces, 4 = Remove non-manifold vertices, 5 = Split non-manifold vertices by threshold.
#' @return cleaned mesh
#' @examples
#' data(humface)
#' cleanface <- humface
#' ##add duplicated faces
#' cleanface$it <- cbind(cleanface$it, cleanface$it[,1:100])
#' ## add duplicated vertices
#' cleanface$vb <- cbind(cleanface$vb,cleanface$vb[,1:100])
#' ## ad unreferenced vertices
#' cleanface$vb <- cbind(cleanface$vb,rbind(matrix(rnorm(18),3,6),1))
#' cleanface <- vcgClean(cleanface, sel=1)
#' @export vcgClean
vcgClean <- function(mesh, sel = 0,tol=0)
    {
        if (!inherits(mesh,"mesh3d"))
            stop("argument 'mesh' needs to be object of class 'mesh3d'")
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        storage.mode(sel) <- "integer"
        tmp <- .Call("Rclean", vb, it, sel, tol)
        tmp$vb <- rbind(tmp$vb,1)
        tmp$normals <- rbind(tmp$normals,1)
        class(tmp) <- "mesh3d"
        return(tmp)
        
    }
