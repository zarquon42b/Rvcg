\name{vcgIntersect}
\alias{vcgIntersect}
\title{check if a mesh is intersected by a set of rays}
\usage{
  vcgIntersect(x, mesh)
}
\arguments{
  \item{x}{a triangular mesh of class 'mesh3d' or a list
  containing vertices and vertex normals (fitting the
  naming of 'mesh3d'.}

  \item{mesh}{triangular mesh to be intersected.}
}
\value{
  list with following items: \item{vb }{3 x n matrix
  containing intersection points} \item{normals }{4 x n
  matrix containing homogenous coordinates of normals at
  intersection points} \item{quality }{integer vector: 1
  indicates that a ray has intersected 'mesh' , while 0
  means not} \item{distance }{numeric vector: distances to
  intersection}
}
\description{
  check if a mesh is intersected by a set of rays
}

