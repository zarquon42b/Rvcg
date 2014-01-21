#include "typedef.h"
#include "RvcgIO.h"
#include <Rcpp.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP RupdateNormals(SEXP _vb, SEXP _it, SEXP _type)
{
  // declare Mesh and helper variables
  int select = Rcpp::as<int>(_type);  
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,_vb,_it);
  /*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  
   
  // update normals
  if (select == 0) {
    tri::UpdateNormal<MyMesh>::PerVertex(m);
  } else {
    tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  }
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  Rcpp::NumericMatrix normals(3,m.vn);
    
  //write back
  
  vi=m.vert.begin();
  SimpleTempData<MyMesh::VertContainer,int> indiceout(m.vert);
 
  for (int i=0;  i < m.vn; i++) {
    if( ! vi->IsD() )	{
      normals(0,i) = (*vi).N()[0];
      normals(1,i) = (*vi).N()[1];
      normals(2,i) = (*vi).N()[2];
    }
    ++vi;
  }
    
  return Rcpp::wrap(normals);
}
 

    
