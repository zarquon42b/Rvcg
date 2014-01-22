\name{vcgUpdateNormals}
\alias{vcgUpdateNormals}
\title{update vertex normals of a triangular meshes or point clouds}
\usage{
vcgUpdateNormals(mesh, type = 0, pointcloud = c(10, 0))
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{type}{select the method to compute per-vertex
  normals: 0=area weighted average of surrounding face
  normals; 1 = angle weighted vertex normals.}

  \item{pointcloud}{add options for normal calculation of
  point clouds.}
}
\value{
mesh with updated/created normals
}
\description{
update vertex normals of a triangular meshes or point
clouds
}
\examples{
data(humface)
humface$normals <- NULL # remove normals
humface <- vcgUpdateNormals(humface)
}

