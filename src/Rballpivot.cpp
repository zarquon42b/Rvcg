// Author: Stefan Schlager
// This is based on code from 

#include "typedefTopo.h"
#include <vcg/complex/algorithms/refine_loop.h>
#include <vcg/complex/algorithms/create/ball_pivoting.h>                                                        
#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rballpivoting(SEXP mesh_, SEXP Radius_, SEXP Clustering_, SEXP CreaseThr_, SEXP deleteFaces_) {
  try {
    TopoMyMesh mesh;
    double Radius = as<double>(Radius_);
    double Clustering = as<double>(Clustering_);
    double CreaseThr = as<double>(CreaseThr_);
    bool deleteFaces = as<bool>(deleteFaces_);
    if (deleteFaces) {
      mesh.fn = 0;
      mesh.face.resize(0);
    }
    Rvcg::IOMesh<TopoMyMesh>::mesh3d2Rvcg(mesh,mesh_);
    tri::BallPivoting<TopoMyMesh> pivot(mesh, Radius, Clustering, CreaseThr);           
    pivot.BuildMesh();
 
    return Rvcg::IOMesh<TopoMyMesh>::RvcgToR(mesh);
 
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

