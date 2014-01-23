#include "typedef.h"
#include "RvcgIO.h"
#include <Rcpp.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;


typedef UpdateTopology<MyMesh>::PEdge SimpleEdge;

RcppExport  SEXP RgetEdge(SEXP vb_, SEXP it_, SEXP unique_)
{
  int i;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  bool unique = Rcpp::as<bool>(unique_);  
 
  // allocate and fill mesh
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  //enable ocf
  /*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  
  // create int indices per face and per vertex to return to R
  SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
  SimpleTempData<MyMesh::FaceContainer,int> indicesf(m.face);
  vi = m.vert.begin();
  for (i=0; i < m.vn; i++) {
    indices[vi] = i;
    ++vi;
  }
  fi = m.face.begin();
  for (i=0; i < m.fn ; i++) {
    indicesf[fi] = i;
    ++fi;
  }
    
  std::vector<SimpleEdge> Edges;
  std::vector< SimpleEdge >::iterator ei;
  std::vector< SimpleEdge >::size_type size;
  tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
  tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m); 
  
  if (unique)
    tri::UpdateTopology<MyMesh>::FillUniqueEdgeVector(m,Edges,true);
  else
    tri::UpdateTopology<MyMesh>::FillEdgeVector(m,Edges,true);
  
  size=Edges.size();
  Rcpp::IntegerVector facept(size), border(size);
  Rcpp::IntegerMatrix edges(size,2);
  border = border * 0;
  
  VertexPointer vp , vp1;
  FacePointer fp;
  i=0;
  
  // for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
  for (i = 0;i < size;i++) {
    vp=Edges[i].v[0];
    vp1=Edges[i].v[1];
    fp=Edges[i].f;
    if( (*fp).IsB(Edges[i].z)==true )
      border[i] = 1;
    
    edges(i,0)=indices[vp]+1;
    edges(i,1)=indices[vp1]+1;
    facept[i] = indicesf[Edges[i].f]+1;
  }
  
  return Rcpp::List::create(Rcpp::Named("edges") = edges,
			    Rcpp::Named("facept") = facept,    
			    Rcpp::Named("border") = border
			    );
}
    
