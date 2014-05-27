#include "typedefCurvature.h"
#include "RvcgIO.h"
//#include <wrap/ply/plylib.cpp>
#include <Rcpp.h>

//using namespace vcg;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP Rcurvature( SEXP vb_, SEXP it_)
{
  // declare Mesh and helper variables
  int i, j;
  CurveMyMesh m;
  VertexIterator vi;
  FaceIterator fi;
 
  Rvcg::IOMesh<CurveMyMesh>::RvcgReadR(m,vb_,it_);
  /*m.vert.EnableQuality();
  m.vert.EnableCurvatureDir();
  m.vert.EnableCurvature();
  m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  
  tri::UpdateTopology<CurveMyMesh>::FaceFace(m);
  tri::UpdateTopology<CurveMyMesh>::VertexFace(m);
  tri::UpdateBounding<CurveMyMesh>::Box(m);
  tri::Allocator<CurveMyMesh>::CompactVertexVector(m);
  tri::UpdateCurvature<CurveMyMesh>::MeanAndGaussian(m);
  tri::UpdateQuality<CurveMyMesh>::VertexFromRMSCurvature(m);
   
  //Bordersearch
  tri::UpdateFlags<CurveMyMesh>::FaceBorderFromNone(m);
  tri::UpdateSelection<CurveMyMesh>::FaceFromBorderFlag(m);
  tri::UpdateFlags<CurveMyMesh>::VertexBorderFromNone(m);
  tri::UpdateSelection<CurveMyMesh>::VertexFromBorderFlag(m);
  
  std::vector<float> gaussvb, meanvb, gaussitmax, meanitmax;
  std::vector<float> RMSvb;
  std::vector<int> bordervb, borderit;
  vi=m.vert.begin();
  //for(i=0; i < m.vn; i++)
  for(i=0; i < m.vn; i++) {
    gaussvb.push_back(vi->Kg());
    meanvb.push_back(vi->Kh());
    RMSvb.push_back(vi->Q());
    if ((*vi).IsS())
      bordervb.push_back(1);
    else
      bordervb.push_back(0);
    ++vi;    
  }
  
  fi=m.face.begin();
  float tmpg, tmpm;
  for(i=0; i < m.fn; i++)  {// get max curvature of vertices per face
    tmpg = (*fi).V(0)->Kg();
    tmpm = (*fi).V(0)->Kh();
    
    for (j = 1; j < 3; j++) {
      if (abs(tmpg) < (*fi).V(j)->Kg())
	tmpg = (*fi).V(j)->Kg();
      if (abs(tmpm) < (*fi).V(j)->Kh())
	tmpm = (*fi).V(j)->Kh();
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
			    Rcpp::Named("meanitmax") = meanitmax
			    );
}
