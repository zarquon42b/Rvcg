#' Import ascii or binary PLY files.
#' 
#' Reads Polygon File Format (PLY) files and stores the results in an object of
#' class "mesh3d" - momentarily only triangular meshes are supported.
#' 
#' 
#' @param file character: file to be read.
#' @param updateNormals logical: if TRUE and the imported file contais faces,
#' vertex normals will be (re)calculated.
#' @param clean logical: if TRUE, duplicated and unreference vertices will be
#' removed.
#' @return Object of class "mesh3d"
#' 
#' with:
#' @return
#' \item{vb }{4xn matrix containing n vertices as homolougous coordinates}
#' \item{normals }{4xn matrix containing vertex normals}
#' \item{it }{4xm matrix containing vertex indices forming triangular faces}
#' \item{material$color }{Per vertex colors if specified in the imported file}
#' @author Stefan Schlager
#' @seealso \code{\link{vcgSmooth}},
#' @keywords ~kwd1 ~kwd2
#' @export 

vcgPlyRead <-function (file,updateNormals=TRUE,clean=TRUE)
{
  ncfile <- nchar(file)
  ext <- substr(file,ncfile-2,ncfile)
  if (ext != "ply" && ext != "PLY")
    stop("please select PLY file")
### get infos from file header
  file <- path.expand(file)
  x <- file
  A <- readLines(x, n = 100, warn=FALSE)
  end <- which(A == "end_header")
  infos <- A[1:end]
  vertinfo <- strsplit(A[grep("element vertex", infos)], " ")
  faceinfo <- strsplit(A[grep("element face", infos)], " ")
  vertbegin <- grep("element vertex",infos)
  facebegin <- grep("element face",infos)
  allvertinfo <- infos[vertbegin:(facebegin-1)]
  allfaceinfo <- infos[facebegin:(end-1)]

  vn <- fn <- 0
  storage.mode(vn) <- "integer"
  storage.mode(fn) <- "integer"
  
  fn <- as.numeric(faceinfo[[1]][3])
  vn <- as.numeric(vertinfo[[1]][3])
  
  
### set flags
  upNorm <- FALSE
  hasQuality <- FALSE
  hasNormals <- FALSE
  hasFaces <- FALSE
  hasVertices <- FALSE
  hasColor <- FALSE
  if (length(grep("property float nx",allvertinfo)) > 0)
    hasNormals <- TRUE
  
  if (fn > 0)
    {
      hasFaces <- TRUE
      if (updateNormals)
        {
          hasNormals <- TRUE
          upNorm <- TRUE
        }
    }
  
  if (vn > 0)
    hasVertices <- TRUE

  if (!hasVertices)
    stop("mesh contains no vertices\n")
  
  color <- grep("property uchar red", infos)
  colvec <- 0
  if (length(color) > 0)
    {
      hasColor <- TRUE
      colvec <- matrix(0,3,vn)
      storage.mode(colvec) <- "integer"
    }
### initialize mesh elements in R
  texinfo<-NULL
  colmat <- NULL
  material <- NULL

  mesh <- list()
  class(mesh) <- "mesh3d"
  vb <- matrix(0,3,vn)
  it <- 0
  quality <- 0
  if (hasFaces)
    it <- matrix(0,3,fn)
  fail <- 0
  storage.mode(it) <- "integer"
  storage.mode(vn) <- "integer"
  storage.mode(fn) <- "integer"
  normals <- 0
  if (hasNormals)
    normals <- vb

### import file ###
  out <- .C("RPlyRead",file,vb,vn,it,fn,normals,as.integer(hasNormals),as.integer(upNorm),quality,as.integer(hasColor),colvec,as.integer(clean),as.integer(fail))

### check if actual vertices or faces exceed numbers read from file header
  if (out[[13]] == 1)
    {
      warning("mesh converted to triangular mesh (maybe Quads involved). Please check result!")
      ## allocate empty mesh with correct numbers 
      vn <- out[[3]]
      fn <- out[[5]]
      vb <- matrix(0,3,vn)
      it <- 0
      quality <- 0
      if (hasFaces)
        it <- matrix(0,3,fn)
      normals <- 0
      if (hasNormals)
        normals <- vb
      out <- .C("RPlyRead",file,vb,vn,it,fn,normals,as.integer(hasNormals),as.integer(upNorm),quality,as.integer(hasColor),colvec,as.integer(clean),as.integer(fail))
    }
  

### fill mesh with imported data
      mesh$vb <- rbind(out[[2]][,1:out[[3]]],1)
      if (hasFaces)
        mesh$it <- out[[4]][,1:out[[5]]]+1
      if (hasNormals)
        mesh$normals <- rbind(out[[6]][,1:out[[3]]],1)
      if (hasColor)
        {
          colvec <- out[[11]][,1:out[[3]]]
          mesh$material <- list()
          colvec <- rgb(colvec[1,],colvec[2,],colvec[3,],maxColorValue=255)
          colfun <- function(x)
            {
              x <- colvec[x]
              return(x)
            }
          mesh$material$color <- matrix(colfun(mesh$it),dim(mesh$it))
        }
    
  return(mesh)
}

