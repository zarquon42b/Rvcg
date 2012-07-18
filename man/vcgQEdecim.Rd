\name{vcgQEdecim}
\alias{vcgQEdecim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Performs Quadric Edge Decimation on triangulare meshes.
}
\description{
Decimates a mesh by adapting the faces of a mesh either to a target face
number, a percentage or an approximate mesh resolution (a.k.a. mean edge
length}
\usage{
vcgQEdecim(mesh, tarface = NULL, percent = NULL, edgeLength = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{Triangular mesh of class "mesh3d"

}
  \item{tarface}{Integer: set number of target faces.
%%     ~~Describe \code{tarface} here~~
}
  \item{percent}{Numeric: between 0 and 1. Set amount of reduction
  relative to existing face number. Overrides tarface argument.
%%     ~~Describe \code{percent} here~~
}
  \item{edgeLength}{Numeric: tries to decimate according to a target
  mean edge length. Under the assumption of regular triangles, the edges
  are half as long by dividing the triangle into 4 regular smaller
  triangles. 
}
}
\details{
This is basically an adaption of the cli tridecimator from vcglib}
\value{
Returns a reduced mesh of class mesh3d.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mesh, tarface = NULL, percent = NULL, edgeLength = NULL) 
{
    vb <- mesh$vb[1:3, ]
    it <- mesh$it - 1
    dimit <- dim(it)[2]
    dimvb <- dim(vb)[2]
    outmesh <- list()
    class(outmesh) <- "mesh3d"
    storage.mode(it) <- "integer"
    if (is.null(tarface) && is.null(percent) && is.null(edgeLength)) 
        stop("please enter decimation option")
    if (!is.null(percent)) {
        if (percent <= 0 || percent > 1) 
            stop("percent must be between 0 and 1")
        tarface <- floor(percent * dimit)
    }
    if (!is.null(edgeLength)) {
        currEL <- vcgMeshres(mesh)
        if (currEL >= edgeLength) 
            warning("edges already shorter than required")
        coef <- (currEL/edgeLength)^2
        print(coef)
        tarface <- floor(coef * dimit)
    }
    storage.mode(tarface) <- "integer"
    tmp <- .C("RQEdecim", vb, ncol(vb), it, ncol(it), tarface, 
        vb)
    outmesh$vb <- rbind(tmp[[1]][, 1:tmp[[2]]], 1)
    outmesh$it <- tmp[[3]][, 1:(tmp[[4]])] + 1
    outmesh$normals <- rbind(tmp[[6]][, 1:(tmp[[2]])], 1)
    if (!is.null(edgeLength)) 
        cat(paste("Mean Edge length is", vcgMeshres(outmesh), 
            "\n"))
    return(outmesh)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
