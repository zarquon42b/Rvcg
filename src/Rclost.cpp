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

void Rclost(double *vb ,int *dim, int *it, int *dimit, double *ioclost, int *clostDim, double *normals, double *dis,int *sign)
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
    const int dref = *clostDim;
    int signo = *sign;
   
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
    VertexPointer ivp[d];
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
    vcg::tri::Allocator<MyMesh>::AddVertices(refmesh,dref);
    VertexPointer ivref[dref];
    vi=refmesh.vert.begin();
     
    for (i=0; i < dref; i++) 
      {
	ivref[i]=&*vi;
	x = ioclost[i*3];
	y = ioclost[i*3+1];
	z = ioclost[i*3+2];
	(*vi).P() = CoordType(x,y,z);
	++vi;
      }
    //--------------------------------------------------------------------------------------//
    //
    //                              INITIALIZE SEARCH STRUCTURES
    //
    // Update the FaceProjection flags needed for projection/distance queries
    // Create a static grid (for fast indexing) and fill it 
    //--------------------------------------------------------------------------------------//
    vcg::tri::Append<MyMesh,MyMesh>::Mesh(outmesh,refmesh);
    tri::UpdateBounding<MyMesh>::Box(m);
    tri::UpdateNormals<MyMesh>::PerFaceNormalized(m);//very important !!!
    tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormals<MyMesh>::NormalizeVertex(m);
    float maxDist = m.bbox.Diag();
    float minDist = 1e-10;
    vcg::tri::FaceTmark<MyMesh> mf; 
    mf.SetMesh( &m );
    vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
    TriMeshGrid static_grid;    
    static_grid.Set(m.face.begin(), m.face.end());
    
    
     for(i=0; i < refmesh.vn; i++)
       {
	 Point3f& currp = refmesh.vert[i].P();
	 Point3f& clost = outmesh.vert[i].P();
	 MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
	 
	 int f_i = vcg::tri::Index(m, f_ptr);
	 MyMesh::CoordType tt = (m.face[f_i].V(0)->N()+m.face[f_i].V(1)->N()+m.face[f_i].V(2)->N());
	 tt=tt/sqrt(tt.dot(tt));
	 dis[i] = minDist;
	 if (signo == 1)
	   {
	     Point3f dif = clost - currp;
	     float sign = dif.dot(tt);	
	     if (sign < 0)
	       { 
		 dis[i] = -dis[i] ;
	       }	
	  }
	 //write back output
	  ioclost[i*3] = clost[0];
	  ioclost[i*3+1] =clost[1];
	  ioclost[i*3+2] =clost[2];
	  normals[i*3] = tt[0];
	  normals[i*3+1] = tt[1];
	  normals[i*3+2] = tt[2];
       }
    //write back output
    /*
      vi=outmesh.vert.begin();
      for (i=0; i< dref; i++) 
      {
	ioclost[i*3] = (*vi).P()[0];
	ioclost[i*3+1] = (*vi).P()[1];
	ioclost[i*3+2] = (*vi).P()[2];
	normals[i*3] = (*vi).N()[0];
	normals[i*3+1] = (*vi).N()[1];
	normals[i*3+2] = (*vi).N()[2];
	++vi;
	}
   
    */
  }
}
