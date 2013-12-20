\name{vcgPlyWrite}
\alias{vcgPlyWrite}
\title{Export meshes to PLY-files}
\usage{
vcgPlyWrite(mesh, filename = dataname, binary = TRUE, addNormals = FALSE,
  writeCol = TRUE)
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{filename}{character: filename (file extension
  '.ply' will be added automatically.}

  \item{binary}{logical: write binary file}

  \item{addNormals}{logical: compute per-vertex normals and
  add to file}

  \item{writeCol}{logical: export existing per-vertex color
  stored in mesh$material$color}
}
\description{
Export meshes to PLY-files (binary or ascii)
}
\examples{
data(humface)
vcgPlyWrite(humface,filename = "humface")
}

