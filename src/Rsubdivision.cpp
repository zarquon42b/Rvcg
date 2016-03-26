// Author: Stefan Schlager
// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include "typedefTopo.h"
#include <vcg/complex/algorithms/refine_loop.h>

#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rsubdivision(SEXP mesh_ ,SEXP iterations_, SEXP threshold_,SEXP type_,SEXP looptype_, SEXP silent_) {
  try { 
    // declare Mesh and helper variables
    int i;
    TopoMyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    float threshold = as<float>(threshold_);
    int iterations = as<int>(iterations_);
    int check = Rvcg::IOMesh<TopoMyMesh>::mesh3d2Rvcg(m,mesh_);
    int type = as<int>(type_);
    int looptype = as<int>(looptype_);
    bool silent = as<bool>(silent_);
    /*m.vert.EnableVFAdjacency();
      m.face.EnableFFAdjacency();
      m.face.EnableVFAdjacency();*/
    if (check != 0) {
      ::Rf_error("mesh has no faces and/or no vertices, nothing done");
    }  else {
      
      tri::UpdateFlags<TopoMyMesh>::FaceBorderFromFF(m);
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      if (tri::Clean<TopoMyMesh>::CountNonManifoldEdgeFF(m) > 0)
	::Rf_error("Mesh has some not 2 manifoldfaces, subdivision surfaces require manifoldness");
      if (threshold < 0) {
	tri::UpdateBounding<TopoMyMesh>::Box(m);
	threshold = m.bbox.Diag()*0.01;
	if (!silent)
	  Rprintf("Edge Threshold set to %f\n",threshold);
      }
      for(int i=0; i< iterations; ++i) {
	if (type == 0) {
	  Refine<TopoMyMesh,MidPointButterfly<TopoMyMesh> >(m,MidPointButterfly<TopoMyMesh>(m),threshold);
	} else if (type == 1) {
	  switch(looptype) {
	  case 0:
	    tri::RefineOddEven<TopoMyMesh>(m, tri::OddPointLoop<TopoMyMesh>(m), tri::EvenPointLoop<TopoMyMesh>(), threshold);
	    break;
	  case 1:
	    tri::RefineOddEven<TopoMyMesh>(m, tri::OddPointLoopGeneric<TopoMyMesh, Centroid<TopoMyMesh>, RegularLoopWeight<TopoMyMesh::ScalarType> >(m),
					   tri::EvenPointLoopGeneric<TopoMyMesh, Centroid<TopoMyMesh>, RegularLoopWeight<TopoMyMesh::ScalarType> >(), threshold);
	    break;
	  case 2:
	    tri::RefineOddEven<TopoMyMesh>(m, tri::OddPointLoopGeneric<TopoMyMesh, Centroid<TopoMyMesh>, ContinuityLoopWeight<TopoMyMesh::ScalarType> >(m),
					   tri::EvenPointLoopGeneric<TopoMyMesh, Centroid<TopoMyMesh>, ContinuityLoopWeight<TopoMyMesh::ScalarType> >(), threshold);
	    break;
        
	  }
	}
      }
    }
  
    return Rvcg::IOMesh<TopoMyMesh>::RvcgToR(m);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}





