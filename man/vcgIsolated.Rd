\name{vcgIsolated}
\alias{vcgIsolated}
\title{Remove isolated pieces from a surface mesh.}
\usage{
vcgIsolated(mesh, facenum = NULL, diameter = NULL)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d".}

  \item{facenum}{integer: all connected pieces with less
  components are removed. If not specified or 0 and
  diameter is NULL, then only the component with the most
  faces is kept.}

  \item{diameter}{numeric: all connected pieces smaller
  diameter are removed removed. 0 removes all component but
  the largest ones. This option overrides 'facenum'.}
}
\value{
returns the reduced mesh.
}
\description{
Remove isolated pieces from a surface mesh, selected by
face number or diameter. Also the option only to keep the
largest piec can be selected
}
\examples{
require(rgl)
data(humface)
cleanface <- vcgIsolated(humface)
}
\author{
Stefan Schlager
}
\seealso{
\code{\link{vcgPlyRead}}
}
\keyword{~kwd1}
\keyword{~kwd2}

