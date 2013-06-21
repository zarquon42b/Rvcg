vcgClean <- function(mesh, sel = 0)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        storage.mode(sel) <- "integer"
        
                                        #tmp <- .C("Rclean",vb=vb,dimvb=ncol(vb),it=it,dimit=ncol(it),norms=vb,sel)
        
        tmp <- .Call("Rclean", vb, it, sel)
        #mesh$vb <- rbind(tmp$vb[,1:tmp$dimvb],1)
        #mesh$it <- tmp$it[,1:tmp$dimit]+1
        #mesh$normals <- tmp$norms[,1:tmp$dimvb]
        tmp$vb <- rbind(tmp$vb,1)
        tmp$normals <- rbind(tmp$normals,1)
        class(tmp) <- "mesh3d"
        return(tmp)
        
    }
