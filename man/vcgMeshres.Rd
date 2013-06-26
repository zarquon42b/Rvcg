\name{vcgMeshres}
\alias{vcgMeshres}
\title{calculate average edge length of a triangular mesh}
\usage{
  vcgMeshres(mesh)
}
\arguments{
  \item{mesh}{triangular mesh stored as object of class
  "mesh3d"}
}
\value{
  Item res average edge length (a.k.a. mesh resolution)

  Item edgelength vector containing lengths for each edge
}
\description{
  calculate average edge length of a triangular mesh, by
  iterating over all faces.
}
\author{
  Stefan Schlager
}
\keyword{~kwd1}
\keyword{~kwd2}

