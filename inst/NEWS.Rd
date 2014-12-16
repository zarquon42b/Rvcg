\name{NEWS}
\title{News for Package Rvcg}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in Rvcg version 0.10 (2014-12-16)}{
  \subsection{New features}{
    \itemize{
      \item added support for reading face and vertex quality from ply files in vcgImport
      \item added option split in vcgIsolated to split a mesh by
  connected components and return them as a list
      \item vcgClean now allows to fix uncoherently oriented faces
      \item added option tol to specify a search radius in vcgClost
      \item added example and docu to vcgUniformRemesh
      \item new function vcgUniformRemesh to remesh based on a voxelized space
      faces for the closest one with a correct normal orientation.
      \item added option weightnorm in vcgClostKD to enable/disable
  calculation of a weighted normal at the closest point
      \item fixes and optimizations in Rkdtree.cpp
      \item added option angdev in vcgClost/Rkdtree to find the closest point with an appropriate normal, if none is found, distance is set to 1e5
      \item added method Laplacian (surface preserving) in vcgSmooth
      \item added vcgStlWrite to export mesh3d objects to STL format.
      \item added option threshold in vcgIsosurface
      \item vcgIsosurface no longer binarizes the array values
    }
  }
}

