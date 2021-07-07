// Author: Stefan Schlager
// Date: 15 September 2010

#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
using namespace Rcpp;

RcppExport SEXP RVFadj(SEXP vb_, SEXP it_)
{
  int i;
  MyMesh m;
  m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();

  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  Rcpp::List outlist(m.vn);
  SimpleTempData<MyMesh::FaceContainer,int>indicesf(m.face);
  typedef vcg::face::VFIterator<MyFace> VFIterator;
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  FaceIterator fi;
  VertexIterator vi;
  fi = m.face.begin();
  for (i = 0; i < m.fn; i++) {
    indicesf[fi] = i;
    ++fi;
  }
  int ptr;
  vi = m.vert.begin();
  for (i = 0; i < m.vn; i++) {
    std::vector<int> tmpvec;
    VFIterator vfi( &*vi);
    while(!vfi.End()) {
      ptr = indicesf[vfi.F()];
      tmpvec.push_back(ptr+1);
      ++vfi;
    }
    outlist[i] = tmpvec;
    ++vi;
  }
  return outlist;
}


// Compute the n-ring vertex neighborhood of the given vertices.
RcppExport SEXP RVVadj(SEXP vb_, SEXP it_, SEXP query_vertices_, SEXP numstep_, SEXP include_self_)
{
  int i;
  int numstep = as<int>(numstep_);
  int include_self = as<int>(include_self_);
  IntegerVector query_vertices(query_vertices_);
  MyMesh m;
  m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();

  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  Rcpp::List outlist(m.vn);
  typedef vcg::face::VFIterator<MyFace> VFIterator;
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  FaceIterator fi;
  VertexIterator vi;

  // Create int vertex indices to return to R.
  SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
  vi = m.vert.begin();
  for (int i=0; i < m.vn; i++) {
    indices[vi] = i;
    ++vi;
  }

  // Start neighborhood computation
  std::vector<std::vector<int>> neighborhoods;
  for(int i=0; i < query_vertices.size(); ++i) {
    int qv = query_vertices[i];
    vi = m.vert.begin()+qv;
    MyVertex* pqv = &*vi;
    std::vector<MyVertex*> neigh;
    vcg::face::VVExtendedStarVF<MyFace>(pqv, numstep, neigh);
    //vcg::face::VVStarVF<MyFace>(pqv, neigh);

    // Collect indices of vertices.
    std::vector<int> neighidx;
    if(include_self) {
      neighidx.push_back(qv + 1); // return 1-based R indices.
    }
    for(int j=0; j<neigh.size(); j++) {
      neighidx.push_back(indices[neigh[j]] + 1); // return 1-based R indices.
    }
    neighborhoods.push_back(neighidx);
  }

  return wrap(neighborhoods);
}




