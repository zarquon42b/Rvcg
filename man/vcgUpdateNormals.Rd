\name{vcgUpdateNormals}
\alias{vcgUpdateNormals}
\title{update vertex normals of a triangular mesh}
\usage{
vcgUpdateNormals(mesh, type = 0)
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{type}{select the method to compute per-vertex
  normals: 0=area weighted average of surrounding face
  normals; 1 = angle weighted vertex normals.}
}
\value{
mesh with updated/created normals
}
\description{
update vertex normals of a triangular mesh
}
\examples{
data(humface)
humface$normals <- NULL # remove normals
humface <- vcgUpdateNormals(humface)
}

