#' Export meshes to PLY-files
#'
#' Export meshes to PLY-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.ply' will be added automatically, if missing.
#' @param binary logical: write binary file
#' @param addNormals logical: compute per-vertex normals and add to file
#' @param writeCol logical: export existing per-vertex color stored in mesh$material$color
#' @param writeNormals write existing normals to file
#' @param \dots additional arguments, currently not used.
#' @examples
#' data(humface)
#' vcgPlyWrite(humface,filename = "humface")
#' ## remove it 
#' unlink("humface.ply")
#' @rdname vcgPlyWrite
#' @export 
vcgPlyWrite <- function(mesh, filename, binary = TRUE, ...) UseMethod("vcgPlyWrite")

#' @rdname vcgPlyWrite
#' @export
vcgPlyWrite.mesh3d <- function(mesh, filename=dataname, binary = TRUE, addNormals = FALSE, writeCol=TRUE, writeNormals=TRUE,...)
{
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    if (!grepl("*\\.ply$",filename,ignore.case=TRUE))
        filename <- paste(filename,".ply",sep="")
    if (!is.null(mesh$material$color) && writeCol==TRUE) {
        ## setup color export
        hasCol <- TRUE
        vn <- ncol(vb)
        col = rep("#FFFFFF", vn)
        if (length(mesh$material$color) != vn) {
            if (!is.null(mesh$it))
                if ((length(mesh$material$color) != length(mesh$it))) {
                    stop("mesh color is not correct")
                } else {
                    tmp1 <- data.frame(it = as.vector(mesh$it))
                    
                    tmp1 <- data.frame(it=1:vn)
                    tmp1$rgb <- as.vector(mesh$material$color)
                    tmp1 <- unique(tmp1)
                    col[tmp1$it] <- tmp1$rgb
                }
        } else 
            col<- mesh$material$color
        colvec <- matrix(col2rgb(col), 3, vn, byrow = F)
        storage.mode(colvec) <- "integer"
    }
    binary <- as.logical(binary)
    addNormals <- as.logical(addNormals)
    mesh$it <- mesh$it-1L
    tmp <- .Call("RMeshWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals,0)
}

#' @rdname vcgPlyWrite
#' @export
vcgPlyWrite.matrix <- function(mesh,filename=dataname, binary = TRUE, addNormals=FALSE, ...) {
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    mm <- list()
    mm$vb <- t(mesh)
    class(mm) <- "mesh3d"
    vcgPlyWrite(mm,filename=filename,binary =binary,writeNormals=FALSE)
}
    
#' Export meshes to STL-files
#'
#' Export meshes to STL-files (binary or ascii)
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.stl' will be added automatically.
#' @param binary logical: write binary file
#' @examples
#' data(humface)
#' vcgStlWrite(humface,filename = "humface")
#' unlink("humface.stl")
#' @rdname vcgStlWrite
#' @export 
vcgStlWrite <- function(mesh, filename=dataname, binary = FALSE) {
   if (!inherits(mesh,"mesh3d"))
        stop("mesh must be of class mesh3d")
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    if (!grepl("*\\.stl$",filename,ignore.case=TRUE))
        filename <- paste(filename,".stl",sep="")
    mesh$it <- mesh$it-1L
    addNormals <- FALSE
    #binary <- FALSE
    writeNormals <- FALSE
    writeCol <- FALSE
    tmp <- .Call("RMeshWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals,3)
}
#' Export meshes to OFF-files
#'
#' Export meshes to OFF-files
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.off' will be added automatically.
#' @examples
#' data(humface)
#' vcgOffWrite(humface,filename = "humface")
#' unlink("humface.off")
#' @rdname vcgOffWrite
#' @export 
vcgOffWrite <- function(mesh, filename=dataname) {
    if (!inherits(mesh,"mesh3d"))
        stop("mesh must be of class mesh3d")
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    if (!grepl("*\\.off$",filename,ignore.case=TRUE))
        filename <- paste(filename,".off",sep="")
    mesh$it <- mesh$it-1L
    addNormals <- FALSE
    binary <- FALSE
    writeNormals <- FALSE
    writeCol <- FALSE
    tmp <- .Call("RMeshWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals,1)
    
}
#' Export meshes to OBJ-files
#'
#' Export meshes to OBJ-files
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.obj' will be added automatically.
#' @param writeNormals write existing normals to file
#' @examples
#' data(humface)
#' vcgObjWrite(humface,filename = "humface")
#' unlink("humface.obj")
#' @export 
vcgObjWrite <- function(mesh, filename=dataname,writeNormals=TRUE) {
    if (!inherits(mesh,"mesh3d"))
        stop("mesh must be of class mesh3d")
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    if (!grepl("*\\.obj$",filename,ignore.case=TRUE))
        filename <- paste(filename,".obj",sep="")
    mesh$it <- mesh$it-1L
    addNormals <- FALSE
    binary <- FALSE
    writeCol <- FALSE
    tmp <- .Call("RMeshWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals,2)
    
}

#' Export meshes to WRL-files
#'
#' Export meshes to WRL-files
#'
#' @param mesh triangular mesh of class 'mesh3d' or a numeric matrix with 3-columns
#' @param filename character: filename (file extension '.wrl' will be added automatically.
#' @param writeCol logical: export existing per-vertex color stored in mesh$material$color
#' @param writeNormals write existing normals to file
#' @examples
#' data(humface)
#' vcgWrlWrite(humface,filename = "humface")
#' unlink("humface.wrl")
#' @export 
vcgWrlWrite <- function(mesh, filename=dataname, writeCol=TRUE,writeNormals=TRUE)
{
    hasCol <- FALSE
    colvec <- matrix(0)
    vb <- mesh$vb[1:3,,drop=FALSE]
    if (!is.matrix(vb))
        stop("mesh has no vertices to write")
    dataname <- deparse(substitute(mesh))
    filename <- path.expand(as.character(filename))
    filename <- paste(filename,".wrl",sep="")
    if (!is.null(mesh$material$color) && writeCol==TRUE) {
        ## setup color export
        hasCol <- TRUE
        vn <- ncol(vb)
        col = rep("#FFFFFF", vn)
        if (!is.null(mesh$it))
            tmp1 <- data.frame(it = as.vector(mesh$it))
        else
            tmp1 <- data.frame(it=1:vn)
        tmp1$rgb <- as.vector(mesh$material$color)
        tmp1 <- unique(tmp1)
        col[tmp1$it] <- tmp1$rgb
        colvec <- matrix(col2rgb(col), 3, vn, byrow = F)
        storage.mode(colvec) <- "integer"
    }

    mesh$it <- mesh$it-1L
    mesh$normals <- mesh$normals
    addNormals <- FALSE
    binary <- FALSE
    tmp <- .Call("RMeshWrite", mesh , binary, addNormals, filename, colvec, hasCol,writeNormals,4)
}
