#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;


 
RcppExport SEXP Rarea(SEXP mesh_, SEXP report_= wrap(true)) {
  try {
    // declare Mesh and helper variables
    int i;
    MyMesh mesh;
    bool report = as<bool>(report_);
    FaceIterator face;
    double area = 0.0;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(mesh,mesh_);
    std::vector<double> faceareas;
    if (report)
      faceareas.resize(mesh.fn);
    int faceind = 0;
    for(face=mesh.face.begin(); face != mesh.face.end(); face++) {
      if(!(*face).IsD()) {
	double tmparea = DoubleArea(*face);
	area += tmparea;
	if (report) {
	  faceareas[faceind] = tmparea/2.0;
	}
	faceind++;
      }
    }
      
    if (!report)
      return wrap(area/2.0);
    else
      return List::create(Rcpp::Named("area") = area,
			  Named("pertriangle") = wrap(faceareas)
			  );
			      
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

