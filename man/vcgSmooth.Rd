\name{vcgSmooth}
\alias{vcgSmooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Smooth a triangular mesh}
\description{
Apply different smoothing algorithms to a triangular mesh.}
\usage{
vcgSmooth(mesh, type = c("taubin", "laplace", "HClaplace", "fujiLaplace", 
    "angWeight"), iteration = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{triangular mesh stored as object of class "mesh3d".
%%     ~~Describe \code{mesh} here~~
}
  \item{type}{character: select smoothing algorithm. Available are
    "taubin", "laplace", "HClaplace", "fujiLaplace", 
    "angWeight" (and any sensible abbreviations).
}
  \item{iteration}{integer: number of iterations to run.
}

}
\details{
The algorithms available are Taubin smoothing, Laplacian smoothing and
an improved version of Laplacian smoothing ("HClaplace"). Also available
are Scale dependent laplacian smoothing ("fujiLaplace") and Laplacian angle weighted smoothing ("angWeight")}
\value{ returns an object of class "mesh3d"
  with:

  \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
  \item{normals }{4xn matrix containing vertex normals.}
  
  \item{quality }{vector: containing distances to target.}
  \item{it }{4xm matrix containing vertex indices forming triangular
  faces.}
}
\references{
  Taubin G. 1995. Curve and surface smoothing without shrinkage. In Computer
Vision, 1995. Proceedings., Fifth International Conference on, pages 852 -
857.

Vollmer J., Mencl R. and Mueller H. 1999. Improved Laplacian Smoothing of
Noisy Surface Meshes. Computer Graphics Forum, 18(3):131-138.

Desbrun, M.; Meyer, M.; Schröder, P. and Barr, A. H. (1999). Implicit fairing of irregular meshes using diffusion and curvature flow,   : 317-324.

}
\author{
Stefan Schlager}
\note{
The additional parameters for taubin smooth are hardcoded to the default
values of meshlab, as they appear to be the least distorting}


\seealso{
 \code{\link{vcgPlyRead}}
}
\examples{

}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
