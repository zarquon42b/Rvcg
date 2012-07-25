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

  
  
extern "C" {

  void Rmeshres(double *vb ,int *dim, int *it, int *dimit, double *res)
  {
    /*typedef typename MyMesh::CoordType CoordType;
    typedef typename MyMesh::ScalarType ScalarType;
    */
    //typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid; 
    typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    ScalarType x,y,z;
    int i;
    MyMesh m;
    MyMesh refmesh;
    MyMesh outmesh;
    // section read from input
    const int d = *dim;
    const int faced = *dimit;
    *res=*res*0;
   
   
    //--------------------------------------------------------------------------------------//
    //
    //                                   PREPROCESS
    // Create meshes,
    // Update the bounding box and initialize max search distance
    // Remove duplicates and update mesh properties
    //--------------------------------------------------------------------------------------//
    vcg::tri::Allocator<MyMesh>::AddVertices(m,d);
    vcg::tri::Allocator<MyMesh>::AddFaces(m,faced);
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
   
    VertexIterator vi=m.vert.begin();
    for (i=0; i < d; i++) 
      {
	ivp[i]=&*vi;
	x = vb[i*3];
	y = vb[i*3+1];
	z=  vb[i*3+2];
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
    Point3f tmp0;
    fi=m.face.begin();
    for (i=0; i < faced ; i++) 
      {
	tmp0 = (*fi).V(0)->P()-(*fi).V(1)->P();
	*res = *res+sqrt(tmp0.dot(tmp0));
	tmp0 = (*fi).V(0)->P()-(*fi).V(2)->P();
	*res = *res+sqrt(tmp0.dot(tmp0));
	tmp0 = (*fi).V(1)->P()-(*fi).V(2)->P();
	*res = *res+sqrt(tmp0.dot(tmp0));
	++fi;
      }
    *res=*res/(faced*3);
  }
}
