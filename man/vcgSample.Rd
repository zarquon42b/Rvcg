\name{vcgSample}
\alias{vcgSample}
\title{Subsample mesh surface}
\usage{
vcgSample(mesh, SampleNum = 100, type = c("km", "pd", "mc"), MCsamp = 20,
  geodes = TRUE, strict = FALSE)
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{SampleNum}{integer Number of sampled points}

  \item{type}{seclect sampling type ("mc"=MonteCarlo
  Sampling, "pd"=PoissonDisk Sampling,"km"=kmean
  clustering)}

  \item{MCsamp}{MonteCarlo sample iterations used in
  PoissonDisk sampling.}

  \item{geodes}{maximise geodesic distance between sample
  points (only for Poisson Disk sampling)}

  \item{strict}{if \code{type="pd"} and the amount of
  coordinates exceeds \code{SampleNum}, the resulting
  coordinates will be subsampled again by kmean clustering
  to reach the requested number.}
}
\value{
sampled points
}
\description{
subsample surface of a triangular mesh
}
\details{
Poisson disk subsampling will not generate the exact amount
of coordinates specified in \code{SampleNum}, depending on
\code{MCsamp} the result wil bee more or less coordinates.
}
\examples{
data(humface)
ss <- vcgSample(humface,SampleNum = 500, type="pd")
\dontrun{
require(rgl)
points3d(ss)
}
}

