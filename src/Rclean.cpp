#include "typedefTopo.h"
#include "RvcgIO.h"
#include <Rcpp.h>

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
  TopoMyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<TopoMyMesh>::RvcgReadR(m,vb_,it_);
  /*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  bool silent = as<bool>(silent_);
  // General cleaning and update of topology
  //tri::UpdateFlags<TopoMyMesh>::VertexBorderFromNone(m);
  //tri::UpdateSelection<TopoMyMesh>::VertexFromBorderFlag(m);
  /*
  tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
  tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
  vcg::tri::UpdateFlags<TopoMyMesh>::FaceBorderFromFF(m);
  vcg::tri::UpdateFlags<TopoMyMesh>::VertexBorderFromFace(m);
  */
  //tri::UpdateFlags<TopoMyMesh>::FaceBorderFromNone(m); 
   
  // do all the cleaning
  for (int i=0; i < select.size();i++) {
    int cnt = select[i];
    if (cnt == 0) { 
      int dupvb = tri::Clean<TopoMyMesh>::RemoveDuplicateVertex(m);
      int dupit = tri::Clean<TopoMyMesh>::RemoveDuplicateFace(m);int unref = tri::Clean<TopoMyMesh>::RemoveUnreferencedVertex(m);
      if (!silent)
	Rprintf("removed %d duplicate faces and %d duplicate vertices\n",dupit,dupvb);
    } else if (cnt == 1) { 
      int unref = tri::Clean<TopoMyMesh>::RemoveUnreferencedVertex(m);
      if (!silent)
	Rprintf("removed %d unreferenced vertices\n",unref);
    } else if (cnt == 2) { 
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::VertexBorderFromFace(m);
      rem = tri::Clean<TopoMyMesh>::RemoveNonManifoldFace(m);
      if (!silent)
	Rprintf("removed %d Non-manifold faces\n",rem);
    } else if (cnt == 3) { 
      rem = tri::Clean<TopoMyMesh>::RemoveDegenerateFace(m);
      if (!silent)
	Rprintf("removed %d degenerate faces\n",rem);
    } else if (cnt == 4) {
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::VertexBorderFromFace(m);
      rem = tri::Clean<TopoMyMesh>::RemoveNonManifoldVertex(m);
      if (!silent)
	Rprintf("removed %d Non-manifold vertices\n",rem);
    } else if (cnt == 5) { 
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::FaceBorderFromFF(m);
      vcg::tri::UpdateFlags<TopoMyMesh>::VertexBorderFromFace(m);
      int split =tri::Clean<TopoMyMesh>::SplitNonManifoldVertex(m,tol);
    if (!silent)
      Rprintf("split %d non-manifold vertices\n",split);
    } else if (cnt == 6) { 
      int merge =tri::Clean<TopoMyMesh>::MergeCloseVertex(m,tol);
      if (!silent)
	Rprintf("merged %d close vertices\n",merge);
    
    } else if (cnt == 7) { 
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
      bool a = false;
      bool b = false;
      tri::Clean<TopoMyMesh>::OrientCoherentlyMesh(m, a ,b);
      //if (!silent)
	//Rprintf("merged %d close vertices\n",merge);
    }
  }
    
  // get a vector of which vertices were removed
  std::vector<int> remvert(m.vert.size());
  std::fill(remvert.begin(), remvert.end(),0);
  vi=m.vert.begin();
  
  int j = 0;
  for (i=0;  i < m.vert.size(); i++) {
    if( vi->IsD() )	{
      remvert[i] = 1;
    }
    ++vi;
  }
  //write back
  vcg::tri::Allocator< TopoMyMesh >::CompactVertexVector(m);
  vcg::tri::Allocator< TopoMyMesh >::CompactFaceVector(m);
  tri::UpdateNormal<TopoMyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<TopoMyMesh>::NormalizePerVertex(m);
  SimpleTempData<TopoMyMesh::VertContainer,int> indiceout(m.vert);
  Rcpp::NumericMatrix vbout(3,m.vn), normals(3,m.vn);
  Rcpp::IntegerMatrix itout(3,m.fn);
  vi=m.vert.begin();
  for (i=0;  i < m.vn; i++) {
    indiceout[vi] = i;
    vbout(0,i) = (*vi).P()[0];
    vbout(1,i) = (*vi).P()[1];
    vbout(2,i) = (*vi).P()[2];
    normals(0,i) = (*vi).N()[0];
    normals(1,i) = (*vi).N()[1];
    normals(2,i) = (*vi).N()[2];
    ++vi;
  }
  FacePointer fp;
  
  fi=m.face.begin();
  j = 0;
  for (i=0; i < m.fn; i++) {
    fp=&(*fi);
    itout(0,i) = indiceout[fp->cV(0)]+1;
    itout(1,i) = indiceout[fp->cV(1)]+1;
    itout(2,i) = indiceout[fp->cV(2)]+1;
    ++fi;
  }
  
  return Rcpp::List::create(Rcpp::Named("vb") = vbout,
			    Rcpp::Named("it") = itout,
			    Rcpp::Named("normals") = normals,
			    Rcpp::Named("remvert") = remvert
			    );
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
 

    
