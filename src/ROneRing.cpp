#include <vcg/simplex/face/pos.h>
#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP ROneRing(SEXP vb_, SEXP it_, SEXP both_ = wrap(false))
{
  try {
    bool both = as<bool>(both_);
    MyMesh m;
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    tri::UpdateTopology<MyMesh>::FaceFace(m);
    tri::UpdateTopology<MyMesh>::VertexFace(m);
    NumericVector area(m.vn);
    NumericVector areaf(m.fn);
    VertexIterator vi = m.vert.begin();
    
    for (int i=0;i < m.vn;i++) {
      //TopoMyVertex v = *vi;
      vcg::face::VFIterator<MyFace> vfi(&*vi); //initialize the iterator tohe first face
      ScalarType test = 0;
      for(;!vfi.End();++vfi) {
	  
	  MyFace* f = vfi.F();
	  test +=  DoubleArea<MyFace>(*f)/2;
	  
	  }
      area[i] = test;
      ++vi;
     
    }
    if (both) {
    FaceIterator fi = m.face.begin();
    for (int i=0;i < m.fn;i++) {
      ScalarType test = 0;
      MyFace* start = (&*fi);
      for (int j = 0; j < 3 ; j++) {
      MyVertex * v = (&*fi)->V(j);
      vcg::face::JumpingPos<MyFace> p((&*fi),j,v);// constructor that takes face, edge and vertex
      do{
	p.NextFE();
	if (!(&*p.F())->IsV()){
	  (&*p.F())->SetV();
	  test += DoubleArea<MyFace>(*p.F())/2;
	  
	}
      } while(p.f!=start);
      }
      areaf[i] = test;
      vcg::tri::UpdateFlags<MyMesh>::FaceClearV(m);

      ++fi;
    }
    }
    return Rcpp::List::create(Rcpp::Named("areaverts") = area,
			      Rcpp::Named("areafaces") = areaf
			      );
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
