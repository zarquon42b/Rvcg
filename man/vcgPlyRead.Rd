\name{vcgPlyRead}
\alias{vcgPlyRead}
\title{
Import ascii or binary PLY files.}
\description{
Reads Polygon File Format (PLY) files and stores the results in an
object of class "mesh3d" - momentarily only triangular meshes are supported.}
\usage{
vcgPlyRead(file, updateNormals = TRUE, clean = TRUE)
}
\arguments{
  \item{file}{character: file to be read.
}
  \item{updateNormals}{logical: if TRUE and the imported file contais
    faces, vertex normals will be (re)calculated.
  }
   \item{clean}{logical: if TRUE, duplicated and unreference vertices will
     be removed.
}
}
\value{
  Object of class "mesh3d"

  with:

  \item{vb }{4xn matrix containing n vertices as homolougous coordinates.}
  \item{normals }{4xn matrix containing vertex normals.}
  \item{it }{4xm matrix containing vertex indices forming triangular faces.}
  
}

\author{
Stefan Schlager}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{vcgSmooth}},
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
