#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include<vcg/complex/algorithms/geodesic.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;
  
  
RcppExport SEXP Rdijkstra(SEXP vb_, SEXP it_, SEXP verts_, SEXP tol_)
{
  try {
  // declare Mesh and helper variables
  //int select = Rcpp::as<int>(type_);  
  IntegerVector verts(verts_);
  int n = verts.length();
  double tol = Rcpp::as<double>(tol_);  
  int i, rem;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  m.vert.EnableVFAdjacency();
  m.vert.EnableQuality();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  std::vector<MyVertex*> seedVec;
  for (int i=0; i < n; i++) {
    vi = m.vert.begin()+verts[i];
    seedVec.push_back(&*vi);
  }
  tri::EuclideanDistance<MyMesh> ed;
  //tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
   //tri::Allocator<MyMesh>::CompactEveryVector(m);
  
   //tri::Geodesic<MyMesh>::Compute(m,seedVec,ed);
   tri::Geodesic<MyMesh>::PerVertexDijkstraCompute(m,seedVec,ed);
   std::vector<float> geodist;
   vi=m.vert.begin();
   for (int i=0; i < m.vn; i++) {
     geodist.push_back(vi->Q());
     ++vi;
   }
  //write back
   
     return wrap(geodist);
  // vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
  // vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
  // tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  // tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  // List out = Rvcg::IOMesh<MyMesh>::RvcgToR(m,true);
  // out.attr("class") = "mesh3d";
  // return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
 

    
