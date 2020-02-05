#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP Rclean(SEXP vb_, SEXP it_, SEXP type_, SEXP tol_, SEXP silent_)
{
  try {
  // declare Mesh and helper variables
  //int select = Rcpp::as<int>(type_);  
  IntegerVector select(type_);
  double tol = Rcpp::as<double>(tol_);  
  int i, rem;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();
  bool silent = as<bool>(silent_);
  // General cleaning and update of topology
  //tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
  //tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
  /*
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(m);
  */
  //tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m); 
   
  // do all the cleaning
  for (int i=0; i < select.size();i++) {
    int cnt = select[i];
    if (cnt == 0) { 
      int dupvb = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
      int dupit = tri::Clean<MyMesh>::RemoveDuplicateFace(m);
      if (!silent)
	Rprintf("removed %d duplicate faces and %d duplicate vertices\n",dupit,dupvb);
    } else if (cnt == 1) { 
      int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
      if (!silent)
	Rprintf("removed %d unreferenced vertices\n",unref);
    } else if (cnt == 2) { 
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      tri::UpdateTopology<MyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(m);
      rem = tri::Clean<MyMesh>::RemoveNonManifoldFace(m);
      if (!silent)
	Rprintf("removed %d Non-manifold faces\n",rem);
    } else if (cnt == 3) { 
      rem = tri::Clean<MyMesh>::RemoveDegenerateFace(m);
      if (!silent)
	Rprintf("removed %d degenerate faces\n",rem);
    } else if (cnt == 4) {
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      tri::UpdateTopology<MyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(m);
      rem = tri::Clean<MyMesh>::RemoveNonManifoldVertex(m);
      if (!silent)
	Rprintf("removed %d Non-manifold vertices\n",rem);
    } else if (cnt == 5) { 
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      tri::UpdateTopology<MyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFaceBorder(m);
      int split =tri::Clean<MyMesh>::SplitNonManifoldVertex(m,tol);
    if (!silent)
      Rprintf("split %d non-manifold vertices\n",split);
    } else if (cnt == 6) { 
      int merge =tri::Clean<MyMesh>::MergeCloseVertex(m,tol);
      if (!silent)
	Rprintf("merged %d close vertices\n",merge);
    
    } else if (cnt == 7) { 
      tri::UpdateTopology<MyMesh>::FaceFace(m);
      tri::UpdateTopology<MyMesh>::VertexFace(m);
      bool a = false;
      bool b = false;
      tri::Clean<MyMesh>::OrientCoherentlyMesh(m, a ,b);
      //if (!silent)
	//Rprintf("merged %d close vertices\n",merge);
    }
  }
    
  // get a vector of which vertices were removed
  std::vector<int> remvert(m.vert.size());
  std::fill(remvert.begin(), remvert.end(),0);
  vi=m.vert.begin();
  
  for (unsigned int j=0;  j < m.vert.size(); j++) {
    if( vi->IsD() )	{
      remvert[j] = 1;
    }
    ++vi;
  }
  //write back
  vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
  vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  List out = Rvcg::IOMesh<MyMesh>::RvcgToR(m,true);
  out["remvert"] = remvert;
  out.attr("class") = "mesh3d";
  return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
 

    
