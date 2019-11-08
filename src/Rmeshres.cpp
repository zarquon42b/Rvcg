#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;


typedef UpdateTopology<MyMesh>::PEdge SimpleEdge;
  
RcppExport SEXP Rmeshres(SEXP vb_ , SEXP it_) {
  try {
    // declare Mesh and helper variables
   
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableVFAdjacency();
    m.face.EnableFFAdjacency();
    m.face.EnableVFAdjacency();
    std::vector<SimpleEdge> Edges;
    std::vector< SimpleEdge >::iterator ei;
    std::vector< SimpleEdge >::size_type size;
    tri::UpdateTopology<MyMesh>::FaceFace(m);
    tri::UpdateTopology<MyMesh>::FillUniqueEdgeVector(m,Edges,true);
    size=Edges.size();
    Rcpp::NumericVector edgelength(size);
    double res = 0;
    double tmp1;
    Point3f tmp0;
    VertexPointer vp , vp1;
    for (std::vector< SimpleEdge >::size_type i = 0;i < size;i++) {
      vp=Edges[i].v[0];
      vp1=Edges[i].v[1];
      tmp0 = vp->P()-vp1->P();
      tmp1 = sqrt(tmp0.dot(tmp0));
      res = res + tmp1;
      edgelength[i] = tmp1;
    }
    res = res/size;
    return Rcpp::List::create(Rcpp::Named("res") = res,
			      Rcpp::Named("edgelength") = edgelength
			      );
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

