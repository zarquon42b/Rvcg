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
#include <cmath>
//#include <utility>  
//#include <algorithm>
extern "C" {

  void Rclost(double *vb ,int *dim, int *it, int *dimit, double *ioclost, int *clostDim, double *normals, double *dis,int *sign,int *border)
  {
    /*typedef MyMesh::CoordType CoordType;
      typedef  MyMesh::ScalarType ScalarType;
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
    typedef MyMesh::VertexPointer VertexPointer;
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
    vcg::tri::Allocator<MyMesh>::AddVertices(refmesh,dref);
    
    //VertexPointer ivref[dref];
    vi=refmesh.vert.begin();
     
    for (i=0; i < dref; i++) 
      {
	
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
    //tri::UpdateNormals<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormals<MyMesh>::PerVertexNormalized(m);
    tri::UpdateNormals<MyMesh>::NormalizeVertex(m);
    float maxDist = m.bbox.Diag()*2;
    float minDist = 1e-10;
    vcg::tri::FaceTmark<MyMesh> mf; 
    mf.SetMesh( &m );
    vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
    TriMeshGrid static_grid;    
    static_grid.Set(m.face.begin(), m.face.end());
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(m);
    
    for(i=0; i < refmesh.vn; i++)
      {
	border[i]=0;
	
	Point3f& currp = refmesh.vert[i].P();
	Point3f& clost = outmesh.vert[i].P();
	MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
	if (f_ptr)
	  {
	    if ((*f_ptr).IsS())
	      border[i]=1;
	     
	    int f_i = vcg::tri::Index(m, f_ptr);
	    MyMesh::CoordType tt;
	    /* ////weighting part momentarily not used because flawed ///
	       std::vector<std::pair<float,int> > xdif;
	       std::vector<float> nweigh(3,0);
	       for (int j=0; j <3;j++)
	       {
	       std::pair <float,int> tmp;
	       tmp=make_pair(0,j);
	       xdif.push_back(tmp);
	       }
		
	       float xsum =0;
	       for (int j=0; j <3;j++)
	       {
	       if (&(m.face[f_i].V(j)->P()))
	       {
	       Point3f vdist = m.face[f_i].V(j)->P() - clost;
	       xdif[j].first = sqrt(vdist.dot(vdist));
	       }
	       xsum = xsum+xdif[j].first;
		
	       }
	       for (int j=0; j <3;j++)
	       {
	       nweigh[j] = xdif[j].first/xsum;//contains weights
	       if (nweigh[j] == 0
	       }
	       //std::sort(nweigh.begin(),nweigh.end(),std::greater<float>());
	       //std::sort(xdif.begin(),xdif.end());*/
	     
	    for (int j=0; j <3;j++)
	      {
		if (&(m.face[f_i].V(j)->N()))
		  {
		    tt +=(m.face[f_i].V(j)->N());
		  }
	      }
	     
	    float vl = sqrt(tt.dot(tt));
	    if (vl > 0 && &vl)//check for zero length normals
	      {
		tt=tt/vl;
	      }   
		 	    
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
      }
  }
}
