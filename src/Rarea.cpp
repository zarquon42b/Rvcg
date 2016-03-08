#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;


 
RcppExport SEXP Rarea(SEXP mesh_) {
  try {
    // declare Mesh and helper variables
    int i;
    MyMesh mesh;
    FaceIterator face;
    double area = 0.0;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(mesh,mesh_);
    for(face=mesh.face.begin(); face != mesh.face.end(); face++)
      if(!(*face).IsD())
	area += DoubleArea(*face);

    return wrap(area/2.0);
			      
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

