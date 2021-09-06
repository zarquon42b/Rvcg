#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include<vcg/complex/algorithms/geodesic.h>

using namespace tri;
using namespace Rcpp;

// Compute pseudo-geodesic distance from query verts_ to all others (or to those
// within a maximal distance of maxdist_ if it is > 0).
RcppExport SEXP Rdijkstra(SEXP vb_, SEXP it_, SEXP verts_, SEXP maxdist_)
{
  try {
    // Declare Mesh and helper variables
    IntegerVector verts(verts_);
    float maxdist = Rcpp::as<float>(maxdist_);
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
    if(maxdist > 0.0) {
      tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed,maxdist);
    } else {
      tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed);
    }

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


RcppExport SEXP RGeodesicPath(SEXP vb_, SEXP it_, SEXP source_, SEXP targets_, SEXP maxdist_)
{
  try {
    // Declare Mesh and helper variables
    int source = Rcpp::as<int>(source_);
    IntegerVector targets(targets_);
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    double maxdist = as<double>(maxdist_);

    // Allocate mesh and fill it
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.vert.EnableQuality();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    tri::UpdateTopology<MyMesh>::VertexFace(m);

    // Create int vertex indices to return to R.
    SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
    vi = m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      indices[vi] = i;
      ++vi;
    }

    // Prepare seed vector with a single vertex
    std::vector<MyVertex*> seedVec;
    vi = m.vert.begin()+source;
    seedVec.push_back(&*vi);

    std::vector<MyVertex*> inInterval;
    MyMesh::PerVertexAttributeHandle<VertexPointer> sourcesHandle;
    sourcesHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m, "sources");
    MyMesh::PerVertexAttributeHandle<VertexPointer> parentHandle;
    parentHandle =  tri::Allocator<MyMesh>::AddPerVertexAttribute<MyMesh::VertexPointer> (m, "parent");

    // Compute pseudo-geodesic distance by summing dists along shortest path in graph.
    tri::EuclideanDistance<MyMesh> ed;
    tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed, maxdist, &inInterval, &sourcesHandle, &parentHandle);
    std::vector<float> geodist;
    vi=m.vert.begin();
    for (int i=0; i < m.vn; i++) {
      geodist.push_back(vi->Q());
      ++vi;
    }

    std::vector<std::vector<int>> paths;
    for(int i=0; i<targets.size(); ++i) {
      int target_vertex = targets[i];
      int current_vertex = target_vertex;
      std::vector<int> path;
      path.push_back(current_vertex + 1);
      while(current_vertex != source) {
        MyMesh::VertexPointer parent = parentHandle[current_vertex];
        int next_vertex = indices[parent];
        current_vertex = next_vertex;
        path.push_back(current_vertex + 1);
      }

      std::reverse(path.begin(), path.end());
      paths.push_back(path);
    }

    List L = List::create(Named("paths") = paths , _["geodist"] = geodist);
    return L;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

