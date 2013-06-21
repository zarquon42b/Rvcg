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
\item{topo }{logical: if TRUE, mesh topology is preserved.
}
\item{quality }{logical: if TRUE, vertex quality is considered.
}
\item{bound }{logical: if TRUE, mesh boundary is preserved.
}
}
\details{
This is basically an adaption of the cli tridecimator from vcglib}
\value{
Returns a reduced mesh of class mesh3d.
}

\author{
Stefan Schlager}


\seealso{
  \code{\link{vcgSmooth}}
  }
\examples{
require(rgl)
data(humface)
##reduce faces to 50\% 
decimface <- vcgQEdecim(humface, percent = 0.5)
## view
shade3d(decimface, col=3)
## some light smoothing
decimface <- vcgSmooth(decimface,iteration = 1)


}


\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
