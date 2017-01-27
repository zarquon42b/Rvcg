// Author: Stefan Schlager
// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include "typedefTopo.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
using namespace vcg;
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Risolated(SEXP vb_ , SEXP it_, SEXP diam_, SEXP facenum_,SEXP silent_, SEXP split_) {
  try { 
    // declare Mesh and helper variables
    int i;
    TopoMyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    bool silent = as<bool>(silent_);
    bool split = as<bool>(split_);
    int check = Rvcg::IOMesh<TopoMyMesh>::RvcgReadR(m,vb_,it_);
    /*m.vert.EnableVFAdjacency();
      m.face.EnableFFAdjacency();
      m.face.EnableVFAdjacency();*/
    if (check != 0) {
      ::Rf_error("mesh has no faces and/or no vertices, nothing done");
    }  else {
      
      double diameter = Rcpp::as<double>(diam_);
      int connect = Rcpp::as<int>(facenum_); 
      //tri::Clean<TopoMyMesh>::RemoveDuplicateVertex(m);
      tri::UpdateTopology<TopoMyMesh>::FaceFace(m);
      tri::UpdateTopology<TopoMyMesh>::VertexFace(m);
      std::pair<int,int> delInfo;
      std::vector< std::pair<int,TopoMyMesh::FacePointer> > CCV;
      int TotalCC=tri::Clean<TopoMyMesh>::ConnectedComponents(m, CCV);
      //Rprintf("%i\n",TotalCC);
      std::vector<float> chunks;
      std::vector<int> chunkface;
      //int CCm = tri::Clean<TopoMyMesh>::ConnectedComponents(m);
      tri::ConnectedComponentIterator<TopoMyMesh> ci;
      for(unsigned int i1=0;i1<CCV.size();++i1) {
	Box3f bb;
	std::vector<TopoMyMesh::FacePointer> FPV;
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
      
      if (!split) {
	if (diameter == 0)
	  diameter = *std::max_element(chunks.begin(),chunks.end());
	if (connect < 0) {
	  delInfo= tri::Clean<TopoMyMesh>::RemoveSmallConnectedComponentsDiameter(m,diameter);
	} else {
	  if (connect == 0)
	    connect = *std::max_element(chunkface.begin(),chunkface.end());
	  delInfo = tri::Clean<TopoMyMesh>::RemoveSmallConnectedComponentsSize(m,connect);
	}
	if (!silent)
	Rprintf("Removed %i connected components out of %i\n", delInfo.second, delInfo.first); 
      tri::Clean<TopoMyMesh>::RemoveUnreferencedVertex(m);
      // get a vector of which vertices were removed
      std::vector<int> remvert(m.vert.size());
      std::fill(remvert.begin(), remvert.end(),0);
      vi=m.vert.begin();
      for (i=0;  i < m.vert.size(); i++) {
	if( vi->IsD() )	{
	  remvert[i] = 1;
	}
	++vi;
      }

      vcg::tri::Allocator< TopoMyMesh >::CompactVertexVector(m);
      vcg::tri::Allocator< TopoMyMesh >::CompactFaceVector(m);
  
      tri::UpdateNormal<TopoMyMesh>::PerVertexAngleWeighted(m);
      tri::UpdateNormal<TopoMyMesh>::NormalizePerVertex(m);
      List out = Rvcg::IOMesh<TopoMyMesh>::RvcgToR(m,true);
      out["remvert"] = remvert;
      out.attr("class") = "mesh3d";
      return out;
	
      } else {
	int length = CCV.size();
	std::list<List> out;
	for(unsigned int i=0; i < CCV.size(); ++i ) {
	  TopoMyMesh destMesh;
	  tri::UpdateSelection<TopoMyMesh>::Clear(m);
	  CCV[i].second->SetS();
	  tri::UpdateSelection<TopoMyMesh>::FaceConnectedFF(m);
	  //tri::UpdateSelection<TopoMyMesh>::Clear(m);
	  tri::UpdateSelection<TopoMyMesh>::VertexFromFaceLoose(m,true);
	  
	  tri::Append<TopoMyMesh,TopoMyMesh>::Mesh(destMesh, m, true);
	  unsigned int facenum = destMesh.fn;
	  tri::UpdateBounding<TopoMyMesh>::Box(destMesh);
	  float mybboxdiam = destMesh.bbox.Diag();
	  if ((connect < 0 && mybboxdiam >= diameter) || (connect >= 0 && facenum >= connect))
	    {
	      tri::UpdateNormal<TopoMyMesh>::PerVertexNormalized(destMesh);				// vertex normals
	      out.push_back(Rvcg::IOMesh<TopoMyMesh>::RvcgToR(destMesh,true));
	    }
	}
	return wrap(out);
      }
	
    } 
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}





