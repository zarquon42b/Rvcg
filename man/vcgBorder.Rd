\name{vcgBorder}
\alias{vcgBorder}
\title{
find all border vertices and faces of a triangular mesh}
\description{
Detect faces and vertices at the borders of a mesh and mark them.
}
\usage{
vcgBorder(mesh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d"
%%     ~~Describe \code{mesh} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{bordervb }{logical: vector containing boolean value for each
    vertex, if it is a border vertex.}
   \item{borderit }{logical: vector containing boolean value for each
    face, if it is a border vertex.}

}


\author{
Stefan Schlager}

\seealso{
 \code{\link{vcgPlyRead}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
