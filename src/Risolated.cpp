// Author: Stefan Schlager
// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include "typedef.h"
#include "RvcgIO.h"
#include <Rcpp.h>
using namespace vcg;
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Risolated(SEXP _vb , SEXP _it, SEXP _diam, SEXP _facenum)
{
  // declare Mesh and helper variables
  int i;
  MyMesh m;
  VertexIterator vi;
  FaceIterator fi;
  
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,_vb,_it);
  /*m.vert.EnableVFAdjacency();
  m.face.EnableFFAdjacency();
  m.face.EnableVFAdjacency();*/
  
  double diameter = Rcpp::as<double>(_diam);
  int connect = Rcpp::as<int>(_facenum); 
  //tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
  tri::UpdateTopology<MyMesh>::FaceFace(m);
  tri::UpdateTopology<MyMesh>::VertexFace(m);
  std::pair<int,int> delInfo;
  std::vector< std::pair<int,MyMesh::FacePointer> > CCV;
  int TotalCC=tri::Clean<MyMesh>::ConnectedComponents(m, CCV);
  Rprintf("%i\n",TotalCC);
  std::vector<float> chunks;
  std::vector<int> chunkface;
  //int CCm = tri::Clean<MyMesh>::ConnectedComponents(m);
  tri::ConnectedComponentIterator<MyMesh> ci;
  for(unsigned int i1=0;i1<CCV.size();++i1) {
    Box3f bb;
    std::vector<MyMesh::FacePointer> FPV;
    for(ci.start(m,CCV[i1].second);!ci.completed();++ci) {
      FPV.push_back(*ci);
      bb.Add((*ci)->P(0));
      bb.Add((*ci)->P(1));
      bb.Add((*ci)->P(2));
    } 
    float diag = bb.Diag();
    chunks.push_back(diag);
    chunkface.push_back(CCV[i1].first);
  }
  
  if (diameter == 0)
    diameter = *std::max_element(chunks.begin(),chunks.end());
  
  if (connect < 0) {
    delInfo= tri::Clean<MyMesh>::RemoveSmallConnectedComponentsDiameter(m,diameter);
  } else {
    if (connect == 0)
      connect = *std::max_element(chunkface.begin(),chunkface.end());
    delInfo = tri::Clean<MyMesh>::RemoveSmallConnectedComponentsSize(m,connect);
  }
      
  Rprintf("Removed %i connected components out of %i\n", delInfo.second, delInfo.first); 
  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
  vcg::tri::Allocator< MyMesh >::CompactVertexVector(m);
  vcg::tri::Allocator< MyMesh >::CompactFaceVector(m);
  
  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
  SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
  Rcpp::NumericMatrix vb(3, m.vn), normals(3, m.vn);
  Rcpp::IntegerMatrix itout(3, m.fn);
  vi=m.vert.begin();
  for (i=0;  i < m.vn; i++) {
    indices[vi] = i;//important: updates vertex indices
    vb(0,i) = (*vi).P()[0];
    vb(1,i) = (*vi).P()[1];
    vb(2,i) = (*vi).P()[2];
    normals(0,i) = (*vi).N()[0];
    normals(1,i) = (*vi).N()[1];
    normals(2,i) = (*vi).N()[2];
    ++vi;
  }
  
  FacePointer fp;
  fi=m.face.begin();
  for (i=0; i < m.fn;i++) {
    fp=&(*fi);
    if( ! fp->IsD() ) {
      itout(0,i) = indices[fp->cV(0)]+1;
      itout(1,i) = indices[fp->cV(1)]+1;
      itout(2,i) = indices[fp->cV(2)]+1;
      ++fi;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vb") = vb,
			    Rcpp::Named("normals") = normals,
			    Rcpp::Named("it") = itout
			    );
}




