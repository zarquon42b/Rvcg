vcgSphere <- function(subdivision = 3) {
    out <- .Call("Rplatonic",subdivision)
    return(out)
}
