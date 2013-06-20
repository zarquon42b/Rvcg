#include <string.h>
#include <vector>
using namespace std;
#include <stdio.h>
#include <cstddef>

// VCG headers for triangular mesh processing
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>

#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>
#include <vcg/complex/algorithms/update/curvature.h>
#include <vcg/complex/algorithms/update/normal.h>

#include <limits>
//#include <../typedef.h>
//#include <wrap/ply/plylib.cpp>
#include <Rcpp.h>

using namespace vcg;
using namespace Rcpp;

class CurvFace;
class CurvEdge;
class CurvVertex;
struct CurvUsedTypes : public UsedTypes<Use<CurvVertex>		::AsVertexType,
					Use<CurvEdge>			::AsEdgeType,
					Use<CurvFace>			::AsFaceType>{};
class CurvEdge : public Edge<CurvUsedTypes>{};
class CurvVertex  : public Vertex< CurvUsedTypes, 
				   vertex::Coord3f, 
				   vertex::BitFlags, 
				   vertex::Normal3f, 
				   vertex::Mark,
				   vertex::Color4b, 
				   vertex::VFAdj,
				   vertex::Curvaturef,
				   vertex::CurvatureDirf,
				   vertex::Qualityf> {};
class CurvFace    : public Face  <CurvUsedTypes, 
				  face::VertexRef,
				  face::BitFlags,
				  face::Mark, 
				  face::Normal3f,
				  face::VFAdj,
				  face::FFAdj> {};

class CurvMesh : public tri::TriMesh< vector<CurvVertex>, vector<CurvFace > >{};
typedef CurvMesh::ScalarType ScalarType;
typedef  CurvMesh::VertexIterator VertexIterator;
typedef  CurvMesh::VertexPointer VertexPointer;
typedef  CurvMesh::FaceIterator FaceIterator;
typedef  CurvMesh::FacePointer FacePointer;
typedef  CurvMesh::CoordType CoordType;
typedef  CurvMesh::ScalarType ScalarType;

typedef CurvMesh::ConstVertexIterator ConstVertexIterator;
typedef CurvMesh::ConstFaceIterator   ConstFaceIterator;
  
  


RcppExport SEXP Rcurvature( SEXP _vb, SEXP _dim, SEXP _it, SEXP _dimit, SEXP _type)
{
        
  ScalarType x,y,z;
  //double x,y,z;
  int i, j;
    
  CurvMesh m;
  // section read from input
    
  const int d =  Rcpp::as<int>(_dim);
  const int faced = Rcpp::as<int>(_dimit);
  const int curvetype = Rcpp::as<int>(_type);
  Rcpp::IntegerVector it(_it);
  //Rcpp::NumericVector vb(_vb);
  Rcpp::NumericMatrix vb(_vb);
   
  //--------------------------------------------------------------------------------------//
  //
  //                                   PREPROCESS
  // Create meshes,
  // Update the bounding box and initialize max search distance
  // Remove duplicates and update mesh properties
  //--------------------------------------------------------------------------------------//
  vcg::tri::Allocator<CurvMesh>::AddVertices(m,d);
  vcg::tri::Allocator<CurvMesh>::AddFaces(m,faced);
  typedef CurvMesh::VertexPointer VertexPointer;
  std::vector<VertexPointer> ivp;
  ivp.resize(d);
    
  VertexIterator vi=m.vert.begin();
  for (i=0; i < d; i++) 
    {
      ivp[i]=&*vi;
      x = (float) vb(0,i);
      y = (float) vb(1,i);
      z=  (float) vb(2,i);
      (*vi).P() = CoordType(x, y, z);
      ++vi;
    } 
    
  int itx,ity,itz;
  FaceIterator fi=m.face.begin();
  for (i=0; i < faced ; i++) 
    {
      itx = it[i*3];
      ity = it[i*3+1];
      itz = it[i*3+2];
      (*fi).V(0)=ivp[itx];
      (*fi).V(1)=ivp[ity];
      (*fi).V(2)=ivp[itz];
      ++fi;
    }
  /*
   int dupvb = tri::Clean<CurvMesh>::RemoveDuplicateVertex(m);
   int dupit = tri::Clean<CurvMesh>::RemoveDuplicateFace(m);
  
   printf("removed %d dupl. Vertices and %d dupl. Faces\n",dupvb,dupit);
  */
  //tri::UpdateQuality<CurvMesh>::VertexConstant(m, 1);  
  //tri::UpdateQuality<CurvMesh>::FaceConstant(m, 1);  
   
  //tri::UpdateFlags<CurvMesh>::FaceBorderFromVF(m);
  //tri::UpdateFlags<CurvMesh>::FaceBorderFromFF(m);
  tri::UpdateTopology<CurvMesh>::FaceFace(m);
  tri::UpdateTopology<CurvMesh>::VertexFace(m);
  tri::UpdateBounding<CurvMesh>::Box(m);
  /*
    tri::UpdateNormal<CurvMesh>::PerFaceNormalized(m);
    tri::UpdateNormal<CurvMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<CurvMesh>::NormalizePerVertex(m);
  */

  tri::Allocator<CurvMesh>::CompactVertexVector(m);
  tri::UpdateCurvature<CurvMesh>::MeanAndGaussian(m);
  tri::UpdateQuality<CurvMesh>::VertexFromRMSCurvature(m);
   //Bordersearch
  tri::UpdateFlags<CurvMesh>::FaceBorderFromNone(m);
  tri::UpdateSelection<CurvMesh>::FaceFromBorderFlag(m);
  tri::UpdateFlags<CurvMesh>::VertexBorderFromNone(m);
  tri::UpdateSelection<CurvMesh>::VertexFromBorderFlag(m);
  
  //tri::UpdateQuality<CurvMesh>::VertexClamp(m,-4,3);
  std::vector<float> gaussvb, meanvb, gaussitmax, meanitmax;
  std::vector<float> RMSvb;
  std::vector<int> bordervb, borderit;
  vi=m.vert.begin();
  //for(i=0; i < m.vn; i++)
  for(i=0; i < m.vn; i++)
    {
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
  for(i=0; i < m.fn; i++)
    {// get max curvature of vertices per face
      tmpg = (*fi).V(0)->Kg();
      tmpm = (*fi).V(0)->Kh();

      for (j = 1; j < 3; j++)
	{
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
