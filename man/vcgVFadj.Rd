\name{vcgVFadj}
\alias{vcgVFadj}
\title{find all faces belonging to each vertex in a mesh}
\usage{
vcgVFadj(mesh)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"}
}
\value{
list containing one vector per vertex containgin the
indices of the adjacent faces
}
\description{
find all faces belonging to each vertex in a mesh and
report their indices
}

