#include "typedef.h"
#include "RvcgIO.h"
#include <Rcpp.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP Rclean(SEXP vb_, SEXP it_, SEXP type_, SEXP tol_)
{
  // declare Mesh and helper variables
  int select = Rcpp::as<int>(type_);  
  double tol = Rcpp::as<double>(tol_);  
  int i, rem;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  /*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  
  // General cleaning and update of topology
  //tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
  //tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
  int dupvb = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  int dupit = tri::Clean<MyMesh>::RemoveDuplicateFace(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  vcg::tri::UpdateFlags<MyMesh>::FaceBorderFromFF(m);
  vcg::tri::UpdateFlags<MyMesh>::VertexBorderFromFace(m);
  Rprintf("removed %d duplicate faces and %d duplicate vertices\n",dupit,dupvb);
  //tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m); 
   
  // do all the cleaning
    
  if (select == 1) { 
    int unref = tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
    Rprintf("removed %d unreferenced vertices\n",unref);
  } else if (select == 2) { 
    rem = tri::Clean<MyMesh>::RemoveNonManifoldFace(m);
    Rprintf("removed %d Non-manifold faces\n",rem);
  } else if (select == 3) { 
    rem = tri::Clean<MyMesh>::RemoveDegenerateFace(m);
    Rprintf("removed %d degenerate faces\n",rem);
  } else if (select == 4) {
    rem = tri::Clean<MyMesh>::RemoveNonManifoldVertex(m);
    Rprintf("removed %d Non-manifold vertices\n",rem);
  } else if (select == 5) { 
    int split =tri::Clean<MyMesh>::SplitNonManifoldVertex(m,tol);
    Rprintf("split %d non-manifold vertices\n",split);
  } else if (select == 6) { 
    int merge =tri::Clean<MyMesh>::MergeCloseVertex(m,tol);
    Rprintf("merged %d close vertices\n",merge);
  } else if (select != 0) {
    Rprintf("unknown parameter\n");
  }
  
  // get a vector of which vertices were removed
  std::vector<int> remvert(m.vert.size());
  std::fill(remvert.begin(), remvert.end(),0);
  vi=m.vert.begin();
  SimpleTempData<MyMesh::VertContainer,int> indiceout(m.vert);
  int j = 0;
  for (i=0;  i < m.vert.size(); i++) {
    if( vi->IsD() )	{
      remvert[i] = 1;
    }
    ++vi;
  }
  //write back
  vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
  vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
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
}
 

    
