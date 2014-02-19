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
cleaned mesh with an additional entry \item{remvert}{vector
of length = number of vertices before cleaning. Entries = 1
indicate that this vertex was removed; 0 otherwise.}
}
\description{
Apply several cleaning algorithms to surface meshes
}
\details{
available options are: 0 = only duplicated vertices and
faces are removed (always applied before cleaning).
\itemize{ \item{1 = remove unreferenced vertices} \item{2 =
Remove non-manifold Faces} \item{3 = Remove degenerate
faces} \item{4 = Remove non-manifold vertices} \item{5 =
Split non-manifold vertices by threshold} \item{6 = merge
close vertices (radius=\code{tol})} }
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

