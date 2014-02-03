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
  diameter are removed removed. \code{diameter = 0} removes
  all component but the largest ones. This option overrides
  the option \code{facenum}.}
}
\value{
returns the reduced mesh.
}
\description{
Remove isolated pieces from a surface mesh, selected by a
minimum amount of faces or of a diameter below a given
threshold. Also the option only to keep the largest piece
can be selected
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

