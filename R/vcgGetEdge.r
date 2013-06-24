vcgGetEdge <- function(mesh,unique=TRUE)
    {
        vb <- mesh$vb[1:3,]
        it <- mesh$it - 1
        dimit <- dim(it)[2]
        dimvb <- dim(vb)[2]
        tmp <- .Call("RgetEdge",vb,it,unique)
        edvert <- tmp$edges
        edge <- data.frame(vert1=edvert[,1])
        edge$vert1 <- edvert[,1]
        edge$vert2 <- edvert[,2]
        edge$facept <- tmp$facept
        edge$border <-  tmp$border
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
