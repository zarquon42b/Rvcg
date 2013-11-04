\name{vcgIntersect}
\alias{vcgIntersect}
\title{check if a mesh is intersected by a set of rays}
\usage{
  vcgIntersect(x, mesh, tol = 0)
}
\arguments{
  \item{x}{a triangular mesh of class 'mesh3d' or a list
  containing vertices and vertex normals (fitting the
  naming of 'mesh3d'.}

  \item{mesh}{triangular mesh to be intersected.}

  \item{tol}{minimum distance to target mesh}
}
\value{
  list with following items: \item{vb }{3 x n matrix
  containing intersection points} \item{normals }{4 x n
  matrix containing homogenous coordinates of normals at
  intersection points} \item{quality }{integer vector
  containing a value for each vertex of \code{x}: 1
  indicates that a ray has intersected 'mesh' , while 0
  means not} \item{distance }{numeric vector: distances to
  intersection}
}
\description{
  check if a mesh is intersected by a set of rays
}
\details{
  project a mesh along a set of given rays (stored as
  normals) onto a target and return the hit points as well
  as information if the target mesh was hit at all. If
  nothing is hit along the ray, the original point's value
  will be retrurned. If the point is already on the surface
  distance will be 0 but be considered as non intersecting
  along its normal.
}

