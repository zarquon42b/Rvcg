\name{vcgIsolated}
\alias{vcgIsolated}
\title{Remove isolated pieces from a surface mesh.
}
\description{
Remove isolated pieces from a surface mesh, selected by face number or
diameter. Also the option only to keep the largest piec can be selected
}
\usage{
vcgIsolated(mesh, facenum = NULL, diameter = NULL)
}
\arguments{
  \item{mesh}{triangular mesh of class "mesh3d".
%%     ~~Describe \code{mesh} here~~
}
  \item{facenum}{integer: all connected pieces with less components are
    removed. If not specified or 0 and diameter is NULL, then only the
   component with the most faces is kept.
%%     ~~Describe \code{facenum} here~~
}
  \item{diameter}{numeric: all connected pieces smaller diameter are removed
    removed. 0 removes all component but the largest ones. This option
    overrides 'facenum'.
%%     ~~Describe \code{diameter} here~~
}
}

\value{returns the reduced mesh.
}

}
\author{Stefan Schlager
%%  ~~who you are~~
}


\seealso{
\code{\link{vcgPlyRead}},\code{\link{rmVertex}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
