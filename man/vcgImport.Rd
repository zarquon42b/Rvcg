\name{vcgImport}
\alias{vcgImport}
\title{Import common mesh file formats.}
\usage{
vcgImport(file, updateNormals = TRUE, readcolor = FALSE, clean = TRUE)
}
\arguments{
  \item{file}{character: file to be read.}

  \item{updateNormals}{logical: if TRUE and the imported
  file contais faces, vertex normals will be
  (re)calculated. Otherwise, normals will be a matrix
  containing zeros.}

  \item{readcolor}{if TRUE, vertex colors will be read if
  available, otherwise all vertices will be colored white.}

  \item{clean}{if TRUE, duplicated and unreferenced
  vertices are removed (be careful when importing point
  clouds).}
}
\value{
  Object of class "mesh3d"

  with: \item{vb }{4 x n matrix containing n vertices as
  homolougous coordinates} \item{it }{3 x m matrix
  containing vertex indices forming triangular faces}
  \item{normals }{4 x n matrix containing vertex normals
  (homologous coordinates)}
}
\description{
  Import common mesh file formats and store the results in
  an object of class "mesh3d" - momentarily only triangular
  meshes are supported.
}
\examples{
data(humface)
vcgPlyWrite(humface)
readit <- vcgImport("humface.ply")
}
\author{
  Stefan Schlager
}
\seealso{
  \code{\link{vcgSmooth}}
}
\keyword{~kwd1}
\keyword{~kwd2}

