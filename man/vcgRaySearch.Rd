\name{vcgRaySearch}
\alias{vcgRaySearch}
\title{check if a mesh is intersected by a set of rays}
\usage{
vcgRaySearch(x, mesh, mintol = 0, maxtol = 1e+15, mindist = FALSE)
}
\arguments{
  \item{x}{a triangular mesh of class 'mesh3d' or a list
  containing vertices and vertex normals (fitting the
  naming of 'mesh3d'.}

  \item{mesh}{triangular mesh to be intersected.}

  \item{mintol}{minimum distance to target mesh}

  \item{maxtol}{maximum distance to search along ray}

  \item{mindist}{search both ways (ray and -ray) and select
  closest point.}
}
\value{
list with following items: \item{vb }{4 x n matrix
containing intersection points} \item{normals }{4 x n
matrix containing homogenous coordinates of normals at
intersection points} \item{quality }{integer vector
containing a value for each vertex of \code{x}: 1 indicates
that a ray has intersected 'mesh' , while 0 means not}
\item{distance }{numeric vector: distances to intersection}
}
\description{
check if a mesh is intersected by a set of rays (stored as
normals)
}
\details{
project a mesh along a set of given rays (stored as
normals) onto a target and return the hit points as well as
information if the target mesh was hit at all. If nothing
is hit along the ray(within the given thresholds), the
ordinary closest point's value will be returned and the
corresponding entry in \code{quality} will be zero.
}
\examples{
data(humface)
#get normals of landmarks
lms <- vcgClost(humface.lm, humface)
# offset landmarks along their normals for a negative amount of -5mm
lms$vb[1:3,] <- lms$vb[1:3,]+lms$normals[1:3,]*-5
intersect <- vcgRaySearch(lms, humface)
\dontrun{
require(Morpho)
require(rgl)
spheres3d(vert2points(lms),radius=0.5,col=3)
plotNormals(lms,long=5)
spheres3d(vert2points(intersect),col=2) #plot intersections
wire3d(humface,col="white")#'
}
}

