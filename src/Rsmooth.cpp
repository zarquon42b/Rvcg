// Author: Stefan Schlager
// Date: 15 September 2010

#include "typedef.h"
#include "RvcgIO.h" 
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP Rsmooth(SEXP _vb, SEXP _it, SEXP _iteration, SEXP _method, SEXP _lambda,  SEXP _mu, SEXP _delta)
{
  int i;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  //set up parameters 
  int iter = Rcpp::as<int>(_iteration);
  int method = Rcpp::as<int>(_method);
  float lambda =  Rcpp::as<float>(_lambda);
  float mu =  Rcpp::as<float>(_mu);
  ScalarType delta = Rcpp::as<double>(_delta);
  //allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,_vb,_it);
       
  if (method == 0) {
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    size_t cnt=tri::UpdateSelection<MyMesh>::VertexFromFaceStrict(m);
    tri::Smooth<MyMesh>::VertexCoordTaubin(m, iter, lambda, mu, cnt>0);
  } else if (method == 1) { 
    tri::Smooth<MyMesh>::VertexCoordLaplacian(m, iter);
  } else if (method == 2) { 
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::Smooth<MyMesh>::VertexCoordLaplacianHC(m, iter);
  } else if (method == 3) { 
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateFlags<MyMesh>::FaceClearB(m);
    tri::Smooth<MyMesh>::VertexCoordScaleDependentLaplacian_Fujiwara(m,iter,delta);
  } else if (method == 4) { 
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateFlags<MyMesh>::FaceClearB(m);
    tri::Smooth<MyMesh>::VertexCoordLaplacianAngleWeighted(m,iter,delta);
  }
  vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);
  vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  Rcpp::NumericMatrix vb(3, m.vn);
  Rcpp::NumericMatrix normals(3, m.vn);
  Rcpp::IntegerMatrix itout(3, m.fn);
  //write back output
  SimpleTempData<MyMesh::VertContainer,int>indices(m.vert);

  // write back updated mesh
  vi=m.vert.begin();
  for (i=0; i < m.vn; i++) {
    indices[vi] = i;
    vb(0,i) = (*vi).P()[0];
    vb(1,i) = (*vi).P()[1];
    vb(2,i) = (*vi).P()[2];
    normals(0,i) = (*vi).N()[0];
    normals(1,i) = (*vi).N()[1];
    normals(2,i) = (*vi).N()[2];
    ++vi;
  }
  
  FacePointer fp;
  fi=m.face.begin();
  for (i=0; i < m.fn; i++) {
    fp=&(*fi);
    if( ! fp->IsD() ) {
      itout(0,i) = indices[fp->cV(0)]+1;
      itout(1,i) = indices[fp->cV(1)]+1;
      itout(2,i) = indices[fp->cV(2)]+1;
      ++fi;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vb") = vb,
			    Rcpp::Named("normals") = normals,
			    Rcpp::Named("it") = itout
			    );
}
  
   


