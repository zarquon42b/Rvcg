#include "typedef.h"
#include "RvcgIO.h"
//#include <wrap/ply/plylib.cpp>
#include <RcppArmadillo.h>

//using namespace vcg;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP Rcurvature( SEXP vb_, SEXP it_)
{
  try {
    // declare Mesh and helper variables
    int i, j;
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
 
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    m.vert.EnableQuality();
      m.vert.EnableCurvatureDir();
      m.vert.EnableCurvature();
      m.vert.EnableVFAdjacency();
      m.face.EnableFFAdjacency();
      m.face.EnableVFAdjacency();
  
    tri::UpdateTopology<MyMesh>::FaceFace(m);
    tri::UpdateTopology<MyMesh>::VertexFace(m);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::Allocator<MyMesh>::CompactVertexVector(m);
    tri::UpdateCurvature<MyMesh>::MeanAndGaussian(m);
    tri::UpdateQuality<MyMesh>::VertexRMSCurvatureFromHGAttribute(m);
    tri::UpdateCurvature<MyMesh>::PrincipalDirections(m);
   
    //Bordersearch
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(m);
    tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
  
    std::vector<float> gaussvb, meanvb, gaussitmax, meanitmax, K1, K2;
    std::vector<float> RMSvb;
    std::vector<int> bordervb, borderit;
    vi=m.vert.begin();
    //for(i=0; i < m.vn; i++)
    auto KH = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<ScalarType> (m, std::string("KH"));
    auto KG = vcg::tri::Allocator<MyMesh>::GetPerVertexAttribute<ScalarType> (m, std::string("KG"));
     
    for(i=0; i < m.vn; i++) {
      gaussvb.push_back(KG[*vi]);
      meanvb.push_back(KH[*vi]);
      RMSvb.push_back(vi->Q());
      K1.push_back(vi->K1());
      K2.push_back(vi->K2());
		   
      if ((*vi).IsS())
	bordervb.push_back(1);
      else
	bordervb.push_back(0);
      ++vi;    
    }
  
    fi=m.face.begin();
    float tmpg, tmpm;
    for(i=0; i < m.fn; i++)  {// get max curvature of vertices per face
      //auto KH = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<ScalarType> (m, std::string("KH"));
      //auto KG = vcg::tri::Allocator<MyMesh>::GetPerFaceAttribute<ScalarType> (m, std::string("KG"));
      tmpg = KG[(*fi).V(0)];
      tmpm = KH[(*fi).V(0)];
      for (j = 1; j < 3; j++) {
	if (abs(tmpg) < abs(KG[(*fi).V(j)]))
	  tmpg = KG[(*fi).V(j)];
	if (abs(tmpm) < abs(KH[(*fi).V(j)]))
	  tmpm = KH[(*fi).V(j)];
      }
      //write borderinfo
      if ((*fi).IsS())
	borderit.push_back(1);
      else
	borderit.push_back(0);
      gaussitmax.push_back(tmpg);
      meanitmax.push_back(tmpm);
      ++fi;    
    }
  
    //return(wrap(curvevb));
    return Rcpp::List::create(Rcpp::Named("gaussvb") = gaussvb,
			      Rcpp::Named("meanvb") = meanvb,
			      Rcpp::Named("RMSvb") = RMSvb,
			      Rcpp::Named("gaussitmax") = gaussitmax,
			      Rcpp::Named("borderit") = borderit,
			      Rcpp::Named("bordervb") = bordervb,
			      Rcpp::Named("meanitmax") = meanitmax,
			      Rcpp::Named("K1") = K1,
			      Rcpp::Named("K2") = K2
			      );
  } catch (std::exception& e) {
    forward_exception_to_r( e );
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  } return R_NilValue; // -Wall

}
