\name{vcgMeshres}
\alias{vcgMeshres}
\title{
calculate average edge length of a triangular mesh}
\description{
calculate average edge length of a triangular mesh, by iterating over
all faces.

}
\usage{
vcgMeshres(mesh)
}
\arguments{
  \item{mesh}{triangular mesh stored as object of class "mesh3d"
  }
}

\value{
  \item{res}{average edge length (a.k.a. mesh resolution)
  }
  \item{edgelength }{vector containing lengths for each edge
  }
}

\author{
Stefan Schlager}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
