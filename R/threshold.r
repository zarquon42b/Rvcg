threshold <- function(array,lowerbound=0,upperbound=0) {
    storage.mode(array) <- "integer"
    nlen <- as.integer(length(array))
    out <- .Call("threshold",array,lowerbound,upperbound)
    return(out)
}
