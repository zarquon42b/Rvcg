\name{vcgImport}
\alias{vcgImport}
\title{Import common mesh file formats.}
\usage{
vcgImport(file, updateNormals = TRUE, readcolor = FALSE)
}
\arguments{
  \item{file}{character: file to be read.}

  \item{updateNormals}{logical: if TRUE and the imported
  file contais faces, vertex normals will be
  (re)calculated. Otherwise, normals will be a matrix
  containing zeros.}

  \item{readcolor}{if TRUE, vertex colors will be read if
  available, otherwise all vertices will be colored white.}
}
\value{
Object of class "mesh3d"

with: \item{vb }{4xn matrix containing n vertices as
homolougous coordinates} \item{it }{4xm matrix containing
vertex indices forming triangular faces} \item{normals
}{4xn matrix containing vertex normals}
}
\description{
Import common mesh file formats and store the results in an
object of class "mesh3d" - momentarily only triangular
meshes are supported.
}
\author{
Stefan Schlager
}
\seealso{
\code{\link{vcgSmooth}},
}
\keyword{~kwd1}
\keyword{~kwd2}

