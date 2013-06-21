\name{vcgClost}
\alias{vcgClost}
\title{
Project coordinates onto a target triangular surface mesh.}
\description{
For a set of 3D-coordinates the closest matches on a target surface are
determined and normals at as well as distances to that point are calculated. 
}
\usage{
vcgClost(x, mesh, sign = TRUE, barycentric = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{k x 3 matrix containing 3D-coordinates or object of class "mesh3d".
}
  \item{mesh}{triangular surface mesh stored as object of class "mesh3d".
}
  \item{sign}{logical: if TRUE, signed distances are returned.
  }
  \item{barycentric }{logical: if TRUE, barycentric coordinates of the
    hit points are returned.
}
}
\details{
}
\value{
  returns an object of class "mesh3d"
  with:

  \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
  \item{normals }{4xn matrix containing vertex normals.}
  
  \item{quality }{vector: containing distances to target.}
  \item{it }{4xm matrix containing vertex indices forming triangular
    faces.Only available, when x is a mesh.}
  \item{barycoords }{3xm Matrix containing barycentric coordinates of
  closest points; only available if barycentric=TRUE.}
}
\references{  Baerentzen, Jakob Andreas. & Aanaes, H., 2002. Generating Signed Distance Fields From Triangle Meshes. Informatics and Mathematical Modelling.
}
\author{
Stefan Schlager}
\note{If large part of the reference mesh are far away from the target
surface, calculation can become very slow. For this case, the function
"closemeshKD" from the package "Morpho" is suggested.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{vcgPlyRead}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
