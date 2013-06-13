rcpptest <- function(mesh)
  
  {
    vb <- mesh$vb[1:3,]
    it <- mesh$it - 1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    storage.mode(it) <- "integer"

    
    bordervb <- rep(0,dimvb)
     borderit <- rep(0,dimit)
    storage.mode(bordervb) <- "integer"
     storage.mode(borderit) <- "integer"
    #tmp <- .Call("Rcpptest")
    tmp <- .Call("Rcpptest",vb,ncol(vb),it,ncol(it))#,bordervb,borderit)
   
    #invisible(list(bordervb=as.logical(tmp[[5]]),borderit=as.logical(tmp[[6]])))
    return(tmp)
  }
