/*#include <string.h>
#include <vector>
using namespace std;
#include <stdio.h>
#include <cstddef>

// VCG headers for triangular mesh processing
#include<vcg/simplex/edge/base.h>
#include<vcg/simplex/vertex/base.h>
#include<vcg/simplex/face/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/create/platonic.h>
#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/quality.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
#include <vcg/complex/algorithms/intersection.h>
#include <vcg/space/index/grid_static_ptr.h>
#include <vcg/space/index/spatial_hashing.h>
#include <vcg/complex/algorithms/closest.h>
#include <vcg/complex/algorithms/smooth.h>
#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>

// VCG File Format Importer/Exporter
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <vcg/complex/algorithms/update/color.h>*/
#include <../typedef.h>
//#include <wrap/ply/plylib.cpp>

#include <Rcpp.h>
extern "C" SEXP Rcpptest( SEXP _vb, SEXP _dim, SEXP _it, SEXP _dimit){

 
  {
    /*typedef MyMesh::CoordType CoordType;
    typedef  MyMesh::ScalarType ScalarType;
    */
    //typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid; 
    // typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    ScalarType x,y,z;
    int i;
    
    MyMesh m;
    MyMesh refmesh;
    MyMesh outmesh;
    // section read from input
    //convert R-data to C++ data
    const int d =  Rcpp::as<int>(_dim);
    const int faced = Rcpp::as<int>(_dimit);
    
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
    //Allocate target
    vcg::tri::Allocator<MyMesh>::AddVertices(m,d);
    vcg::tri::Allocator<MyMesh>::AddFaces(m,faced);
    typedef MyMesh::VertexPointer VertexPointer;
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
    
    VertexIterator vi=m.vert.begin();
    for (i=0; i < d; i++) 
      {
	ivp[i]=&*vi;
	/*x = vb[i*3];
	y = vb[i*3+1];
	z=  vb[i*3+2];*/
	x = vb[1,i];
	y = vb[2,i];
	z=  vb[3,i];
	(*vi).P() = CoordType(x,y,z);
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
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(m);
    tri::UpdateFlags<MyMesh>::VertexBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::VertexFromBorderFlag(m);
    //std::vector<std::map<std::string,int> > v;
    
    //int i;
    std::vector<int> bordervb, borderit;
    //Rcpp::NumericVector xa();
    
    vi=m.vert.begin();
    for(i=0; i < m.vn; i++)
      {
	
	//bordervb[i]=0;
	
	if ((*vi).IsS())
	   bordervb.push_back(1);
	else
	  bordervb.push_back(0);
	++vi;    
      }
    
    //return Rcpp::wrap(xa);
    fi=m.face.begin();
    for(i=0; i < m.fn; i++)
      { 
	 
	
	if ((*fi).IsS())
	  borderit.push_back(1);
	else
	  borderit.push_back(0);
	 ++fi;    
      }
    return Rcpp::List::create(Rcpp::Named("bordervb") = bordervb,
			      Rcpp::Named("borderit") = borderit
			      );
  }
}
