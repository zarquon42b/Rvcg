vcgGetEdge <- function(mesh,unique=TRUE)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        storage.mode(it) <- "integer"
        edges <- matrix(0,2,3*ncol(it))
        storage.mode(edges) <- "integer"
        unique <- as.integer(unique)
        edgecount <- 0;
        storage.mode(edgecount) <- "integer"    
        tmp <- .C("RgetEdge",vb,ncol(vb),it,ncol(it),edgecount,edges,unique)
        edges <- tmp[[6]][,1:tmp[[5]]]
        edges <- edges+1
        if (!unique)
            edges <- edges[,order(edges[1,],edges[2,])]
        invisible(edges)
    }
    
