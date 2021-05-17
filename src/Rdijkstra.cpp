#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include<vcg/complex/algorithms/geodesic.h>

using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rdijkstra(SEXP vb_, SEXP it_, SEXP verts_)
{
  try {
    // Declare Mesh and helper variables
    IntegerVector verts(verts_);
    int n = verts.length();
    int i, rem;
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    
    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);
    
    // Prepare seed vector
    std::vector<MyVertex*> seedVec;
    for (int i=0; i < n; i++) {
      vi = m.vert.begin()+verts[i];
      seedVec.push_back(&*vi);
    }
    
    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    tri::EuclideanDistance<MyMesh> ed;
    tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed);
    std::vector<float> geodist;
    vi=m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      geodist.push_back(vi->Q());
      ++vi;
    }
    return wrap(geodist);
    
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
