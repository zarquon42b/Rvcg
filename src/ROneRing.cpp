#include <vcg/simplex/face/pos.h>
#include "typedefTopo.h"
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
    TopoMyMesh m;
    Rvcg::IOMesh<TopoMyMesh>::RvcgReadR(m,vb_,it_);
    tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
    tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
    NumericVector area(m.vn);
    NumericVector areaf(m.fn);
    VertexIterator vi = m.vert.begin();
    
    for (int i=0;i < m.vn;i++) {
      //TopoMyVertex v = *vi;
      vcg::face::VFIterator<TopoMyFace> vfi(&*vi); //initialize the iterator tohe first face
      ScalarType test = 0;
      for(;!vfi.End();++vfi) {
	  
	  TopoMyFace* f = vfi.F();
	  test +=  DoubleArea<TopoMyFace>(*f)/2;
	  
	  }
      area[i] = test;
      ++vi;
     
    }
    if (both) {
    FaceIterator fi = m.face.begin();
    for (int i=0;i < m.fn;i++) {
      ScalarType test = 0;
      TopoMyFace* start = (&*fi);
      for (int j = 0; j < 3 ; j++) {
      TopoMyVertex * v = (&*fi)->V(j);
      vcg::face::JumpingPos<TopoMyFace> p((&*fi),j,v);// constructor that takes face, edge and vertex
      do{
	p.NextFE();
	if (!(&*p.F())->IsV()){
	  (&*p.F())->SetV();
	  test += DoubleArea<TopoMyFace>(*p.F())/2;
	  
	}
      } while(p.f!=start);
      }
      areaf[i] = test;
      vcg::tri::UpdateFlags<TopoMyMesh>::FaceClearV(m);

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
