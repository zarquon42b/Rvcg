% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcgSmooth.r
\name{vcgSmoothImplicit}
\alias{vcgSmoothImplicit}
\title{Implicit Smoothes a triangular mesh}
\usage{
vcgSmoothImplicit(
  mesh,
  lambda = 0.2,
  useMassMatrix = TRUE,
  fixBorder = FALSE,
  useCotWeight = FALSE,
  degree = 1L,
  lapWeight = 1,
  SmoothQ = FALSE
)
}
\arguments{
\item{mesh}{triangular mesh stored as object of class "mesh3d".}

\item{lambda}{numeric: the amount of smoothness, useful only if
\code{useMassMatrix} is \code{TRUE}; default is \code{0.2}}

\item{useMassMatrix}{logical: whether to use mass matrix to keep the mesh
close to its original position (weighted per area distributed on vertices);
default is \code{TRUE}}

\item{fixBorder}{logical: whether to fix the border vertices of the mesh;
default is \code{FALSE}}

\item{useCotWeight}{logical: whether to use cotangent weight; default is
\code{FALSE} (using uniform 'Laplacian')}

\item{degree}{integer: degrees of 'Laplacian'; default is \code{1}}

\item{lapWeight}{numeric: weight when \code{useCotWeight} is \code{FALSE};
default is \code{1.0}}

\item{SmoothQ}{logical: whether to smooth the quality (distances to target).}
}
\value{
returns an object of class "mesh3d" with:
\item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
\item{normals}{4xn matrix containing vertex normals.}
\item{it }{4xm matrix containing vertex indices forming triangular
faces.}
}
\description{
Applies implicit smoothing algorithms on a triangular mesh.
}
\examples{

data(humface)
smoothface <- vcgSmoothImplicit(humface)
## view
\dontrun{
require(rgl)
shade3d(smoothface, col=3)
}

}
\seealso{
\code{\link{vcgPlyRead},\link{vcgClean},\link{vcgSmooth}}
}
\author{
Zhengjia Wang
}
