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

RcppExport SEXP Rmeshvol(SEXP mesh_) {
  try {
    MyMesh m;
    int check = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(m,mesh_);
    bool Watertight, Oriented = false;
    int VManifold, FManifold;
    float Volume = 0;
    int numholes, BEdges = 0;
    //check manifoldness
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    m.face.EnableNormal();
    UpdateTopology<MyMesh>::FaceFace(m);
    VManifold = Clean<MyMesh>::CountNonManifoldVertexFF(m);
    FManifold = Clean<MyMesh>::CountNonManifoldEdgeFF(m);
    
    if ((VManifold>0) || (FManifold>0)) {
      ::Rf_error(
        (
          "Mesh is not manifold\n  Non-manifold vertices: " +
          std::to_string(VManifold) +"\n" +
          "  Non-manifold edges: " + 
          std::to_string(FManifold) +"\n"
        ).c_str()
      );
    }
      
     
    Watertight = Clean<MyMesh>::IsWaterTight(m);
    Oriented = Clean<MyMesh>::IsCoherentlyOrientedMesh(m);
    tri::Inertia<MyMesh> mm(m);
    mm.Compute(m);
    Volume = mm.Mass();

    // the sign of the volume depend on the mesh orientation
    if (Volume < 0.0)
      Volume = -Volume;
    if (!Watertight)
      ::Rf_warning("Mesh is not watertight! USE RESULT WITH CARE!");
    if (!Oriented)
      ::Rf_warning("Mesh is not coherently oriented! USE RESULT WITH CARE!");

  return wrap(Volume);

} catch (std::exception& e) {
  ::Rf_error( e.what());
 } catch (...) {
  ::Rf_error("unknown exception");
 }
}
