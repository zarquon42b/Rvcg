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
  
   


