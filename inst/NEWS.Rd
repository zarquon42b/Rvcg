\name{NEWS}
\title{News for Package Rvcg}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}


\section{Changes in Rvcg version 0.18 (2018-09-28)}{
  \subsection{small changes and fixes}{
    \itemize{
      \item Rintersect.cpp: check for valid pointer and increased search radius to 2*bbox
      \item vcgBox: added default argument mesh=vcgSphere()
      \item vcgArea: added option pertriangle to get per face area
      \item vcglib: fixed potential buffer overflow in src/vcglib/wrap/io_trimesh/export_stl.h
    }
  }
}


\section{Changes in Rvcg version 0.17 (2017-09-04)}{
  \subsection{small changes and fixes}{
    \itemize{
      \item vcgClostKD now sets all distances > 1e12 to NaN
      \item vcgOffWrite/vcgObjWrite: fixed typo in documentation
    }
  }
}


\section{Changes in Rvcg version 0.16 (2017-04-04)}{
  \subsection{new features}{
    \itemize{
      \item updated vcglib to version 1.0.1
      \item added meshInfo, nfaces and nverts to obtain mesh information
      \item vcgStlWrite and vcgPlyWrite: avoid doubled file extension if
      filename already contains extension
      \item NEW: vcgObjWrite, vcgOffWrite and vcgWrlWrite, allowing
      export of meshes into OBJ, OFF and WRL format
      \item vcgImport: now supports OFF file format
    }
  }
}

\section{Changes in Rvcg version 0.15 (2017-01-13)}{
  \subsection{new features}{
    \itemize{
      \item updated vcglib to revision 5760 and adapted patches and vcgQEdecim
      \item added meshInfo, nfaces and nverts to obtain mesh information
    }
  }
}
  

\section{Changes in Rvcg version 0.14 (2016-07-06)}{
  \subsection{new features}{
    \itemize{
      \item added vcgBallPivoting: ball pivoting surface reconstruction
      \item added vcgKmeans a fast kmeans clustering for 2D and 3D point clouds
      \item added vcgClostOnKDtreeFromBarycenters, vcgCreateKDtree,
  vcgCreateKDtreeFromBarycenters, vcgSearchKDtree to allow reusing KD-trees
      \item added vcgSubdivide to refine an existing mesh by subidviding faces
      \item added vcgMetro to compare differences between two meshes
  using different types of subsampling (thanks to F. Girinion)
      \item added vcgArea to calculate the surface area of a mesh
      \item added option keep to vcgIsolated and allowed options diameter/facenum to be used with split=T
      \item updated vcglib to revision 5735
      \item added Rvcg::IOMesh::mesh3d2Rvcg<T>, allowing easier
      conversion from mesh3d to vcglib meshtypes
      \item added checkFaceOrientation to check whether a the outer
  layer of a surface mesh points outward
      }
  }
  
  \subsection{bug fixes}{
    \itemize{
      \item vcgClost: in case distances are beyond threshold, distance values are
      to NaN and a  warning is issued
      \item refactored and simplified C++ code

      
      }
  }
}

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

