vcgGetEdge <- function(mesh)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        edges <- matrix(0,2,2*ncol(it))
        storage.mode(edges) <- "integer"
        edgecount <- 0;
        storage.mode(edgecount) <- "integer"    
        tmp <- .C("RgetEdge",vb,ncol(vb),it,ncol(it),edgecount,edges)
        edges <- tmp[[6]][,1:tmp[[5]]]
        edges <- edges+1
        invisible(edges)
    }
    
