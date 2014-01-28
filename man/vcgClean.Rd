\name{vcgClean}
\alias{vcgClean}
\title{Clean triangular surface meshes}
\usage{
vcgClean(mesh, sel = 0, tol = 0)
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{sel}{integer select cleaning type (see "details")}

  \item{tol}{numeric value determining Vertex Displacement
  Ratio used for splitting non-manifold vertices.}
}
\value{
  cleaned mesh
}
\description{
  Apply several cleaning algorithms to surface meshes
}
\details{
  available options are: 0 = only duplicated vertices and
  faces are removed (always applied before cleaning). 1 =
  remove unreferenced vertices, 2 = Remove non-manifold
  Faces, 3 = Remove degenerate faces, 4 = Remove
  non-manifold vertices, 5 = Split non-manifold vertices by
  threshold.
}
\examples{
data(humface)
cleanface <- humface
##add duplicated faces
cleanface$it <- cbind(cleanface$it, cleanface$it[,1:100])
## add duplicated vertices
cleanface$vb <- cbind(cleanface$vb,cleanface$vb[,1:100])
## ad unreferenced vertices
cleanface$vb <- cbind(cleanface$vb,rbind(matrix(rnorm(18),3,6),1))
cleanface <- vcgClean(cleanface, sel=1)
}

