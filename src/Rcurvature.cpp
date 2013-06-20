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
    int i;
    
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
      
    //tri::UpdateQuality<CurvMesh>::VertexConstant(m, 1);  
    //tri::UpdateQuality<CurvMesh>::FaceConstant(m, 1);  
   
    //tri::UpdateFlags<CurvMesh>::FaceBorderFromVF(m);
    //tri::UpdateFlags<CurvMesh>::FaceBorderFromFF(m);
    tri::UpdateTopology<CurvMesh>::FaceFace(m);
    tri::UpdateTopology<CurvMesh>::VertexFace(m);
    tri::UpdateBounding<CurvMesh>::Box(m);
    tri::UpdateNormal<CurvMesh>::PerFaceNormalized(m);
    tri::UpdateNormal<CurvMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<CurvMesh>::NormalizePerVertex(m);

    tri::Allocator<CurvMesh>::CompactVertexVector(m);
    tri::UpdateCurvature<CurvMesh>::MeanAndGaussian(m);
    //tri::UpdateCurvature<CurvMesh>::ComputeSingleVertexCurvature(m);
    
    // if (curvetype == 1)
    {
      // tri::UpdateQuality<CurvMesh>::VertexFromGaussianCurvatureHG(m);
      //tri::UpdateQuality<CurvMesh>::FaceFromRMSCurvature(m);
    }
    //tri::UpdateQuality<CurvMesh>::VertexClamp(m,-4,3);
    std::vector<float> curvevb;
   float cc;
   vi=m.vert.begin();
   //for(i=0; i < m.vn; i++)
 for(i=0; i < m.vn; i++)
      {
	cc = (*vi).Kg();
	// m.vert[i].Kg();
	curvevb.push_back(cc);
	if (i == 1)
	  printf("%f\n",cc);
	++vi;    
      }
   //printf("%d",
   //return(wrap(curvevb));
 return Rcpp::List::create(Rcpp::Named("curvevb") = curvevb
			      
			      );

}
