// Author: Stefan Schlager
// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include "typedef.h"

#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rsubdivision(SEXP mesh_ ,SEXP iterations_, SEXP threshold_,SEXP type_,SEXP looptype_, SEXP silent_) {
  try { 
    // declare Mesh and helper variables
    int i;
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    float threshold = as<float>(threshold_);
    int iterations = as<int>(iterations_);
    int check = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(m,mesh_);
    int type = as<int>(type_);
    int looptype = as<int>(looptype_);
    bool silent = as<bool>(silent_);
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    if (check != 0) {
      ::Rf_error("mesh has no faces and/or no vertices, nothing done");
    }  else {
      
      tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      if (tri::Clean<MyMesh>::CountNonManifoldEdgeFF(m) > 0)
	::Rf_error("Mesh has some not 2 manifoldfaces, subdivision surfaces require manifoldness");
      if (threshold < 0) {
	tri::UpdateBounding<MyMesh>::Box(m);
	threshold = m.bbox.Diag()*0.01;
	if (!silent)
	  Rprintf("Edge Threshold set to %f\n",threshold);
      }
      for(int i=0; i< iterations; ++i) {
	if (type == 0) {
	  Refine<MyMesh,MidPointButterfly<MyMesh> >(m,MidPointButterfly<MyMesh>(m),threshold);
	} else if (type == 1) {
	  switch(looptype) {
	  case 0:
	    tri::RefineOddEven<MyMesh>(m, tri::OddPointLoop<MyMesh>(m), tri::EvenPointLoop<MyMesh>(), threshold);
	    break;
	  case 1:
	    tri::RefineOddEven<MyMesh>(m, tri::OddPointLoopGeneric<MyMesh, Centroid<MyMesh>, RegularLoopWeight<MyMesh::ScalarType> >(m),
					   tri::EvenPointLoopGeneric<MyMesh, Centroid<MyMesh>, RegularLoopWeight<MyMesh::ScalarType> >(), threshold);
	    break;
	  case 2:
	    tri::RefineOddEven<MyMesh>(m, tri::OddPointLoopGeneric<MyMesh, Centroid<MyMesh>, ContinuityLoopWeight<MyMesh::ScalarType> >(m),
					   tri::EvenPointLoopGeneric<MyMesh, Centroid<MyMesh>, ContinuityLoopWeight<MyMesh::ScalarType> >(), threshold);
	    break;
        
	  }
	}
      }
    }
  
    return Rvcg::IOMesh<MyMesh>::RvcgToR(m);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}





