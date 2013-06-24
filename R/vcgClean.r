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
