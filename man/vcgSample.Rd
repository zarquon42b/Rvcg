\name{vcgSample}
\alias{vcgSample}
\title{Subsample mesh surface}
\usage{
  vcgSample(mesh, SampleNum = 10, type = 1, MCsamp = 20,
    geodes = TRUE)
}
\arguments{
  \item{mesh}{triangular mesh of class 'mesh3d'}

  \item{SampleNum}{integer Number of sampled points}

  \item{type}{seclect sampling type (1=MonteCarlo Sampling,
  2=PoissonDisk Sampling)}

  \item{MCsamp}{MonteCarlo sample iterations used in
  PoissonDisk sampling.}

  \item{geodes}{maximise geodesic distance between sample
  points}
}
\value{
  sampled points
}
\description{
  subsample surface of a triangular mesh
}
\details{
  not ready yet
}
\examples{
require(rgl)
data(humface)
ss <- vcgSample(humface,SampleNum = 500, type=2)
points3d(ss)
}

