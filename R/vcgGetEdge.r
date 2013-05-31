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
        border <- facept <- as.integer(rep(0,3*ncol(it)))
        storage.mode(edgecount) <- "integer"    
        tmp <- .C("RgetEdge",vb,ncol(vb),it,ncol(it),edgecount,edges,facept,border,unique)
        
        edvert <- t(tmp[[6]][,1:tmp[[5]]]+1)
        
        edge <- data.frame(vert1=edvert[,1])
        #edge$vert1 <- edvert[,1]
        edge$vert2 <- edvert[,2]
        edge$facept <- (tmp[[7]][1:tmp[[5]]])+1
        
        edge$border <-  (tmp[[8]][1:tmp[[5]]])
        if (!unique)
            edge <- edge[order(edge[,1],edge[,2]),]
        invisible(edge)
    }
    
vcgNonBorderEdge <- function(mesh)
    {
        edges <- vcgGetEdge(mesh,unique=FALSE)
        border <- which(edges$border == 1)
        edgesClean <- edges
        if (length(border) != 0)
            edgesClean <- edgesClean[-border,]

        cat(paste("mesh contains ",length(border), "border edges\n"))
        n <- dim(edgesClean)[1]/2
        n2 <- (1:n)*2
        n1 <- n2-1
        out <- edgesClean[n1,]
        out$face1 <- edgesClean$facept[n1]
        out$face2 <- edgesClean$facept[n2]

        return(out)
    }
