\name{NEWS}
\title{News for Package Rvcg}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in Rvcg version 0.13.1.2 (2016-01-22)}{
  \subsection{bug fixes}{
    \itemize{
      \item  changed SUPPORT_OPENMP to _OPENMP in C++ code to comply
      with CRAN rules.
      }
  }
}

\section{Changes in Rvcg version 0.13.1.1 (2015-12-31)}{
  \subsection{bug fixes}{
    \itemize{
      \item fixed clang build for OS X and Fedora 
      }
  }
}

\section{Changes in Rvcg version 0.13.1 (2015-12-17)}{
  \subsection{bug fixes}{
    \itemize{
      \item fixed build for solaris
      }
  }
}
  
\section{Changes in Rvcg version 0.13.0 (2015-12-15)}{
  \subsection{new features}{
    \itemize{
      \item updated vcglib source code to Revision: 5542
      \item added functions to create basic meshes (spheres,
      polyhedrons)
      \item enabled OpenMP support to run closest point search in
      parallel (vcgKDtree and vcgClostKD)
      \item added options IJK2RAS and direction in vcgIsosurface
      \item added new option writeNormals in vcgPlyWrite
      \item added option facenormals in vcgClost* functions
      \item added functions to check if the normals of a mesh are
      oriented outward
      \item added additional low-level options in vcgKDtree

    }
  }
  \subsection{bug fixes}{
    \itemize{
      \item fixed useless check in vcgUpdateNormals
      \item added BugReport page in DESCRIPITON
    }
  }
}

\section{Changes in Rvcg version 0.12.0 (2015-06-23)}{
  \subsection{new features}{
    \itemize{
      \item updated vcglib source code to Revision: 5521
      \item added option as.int to vcgIsosurface
      \item added option iterate in vcgClean
    }
  }
  \subsection{bug fixes}{
    \itemize{
      \item removed wrong error call in vcgKDtree
      \item silenced all output in vcgQEdecim if silent=T
    }
  }
}


\section{Changes in Rvcg version 0.11 (2015-04-14)}{
  \subsection{new features}{
    \itemize{
      \item updated vcglib source code to Revision: 5496
      \item added support for xyz files
      \item improved handling of texture files
      }
  }
  \subsection{bug fixes}{
    \itemize{
      \item fixed Rkdtree.cpp where a flying snippet could provoke assertion errors
      \item fixed issues if target of a vcgKdtree query is a mesh
      \item disabled cleaning for pointclouds in RallRead.cpp      
    }
  }
}


\section{Changes in Rvcg version 0.10.1 (2014-12-18)}{
  \subsection{bug fixes}{
    \itemize{
      \item fixed texture import in vcgImport
      \item fixed missing <#include time.h> in upstream code preventing
      windows build
      
    }
  }
}


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

