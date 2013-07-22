#include "typedef.h"
//#include <wrap/ply/plylib.cpp>

  
  
extern "C" {

void Rintersect(double *vb ,int *dim, int *it, int *dimit, double *ioclost, int *clostDim, double *normals, double *dis, int *hitbool)
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
    
    
    float t;
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
    Point3f normtmp;
    for (i=0; i < dref; i++) 
      {
	
	x = ioclost[i*3];
	y = ioclost[i*3+1];
	z = ioclost[i*3+2];
	(*vi).P() = CoordType(x,y,z);
	x = normals[i*3];
	y = normals[i*3+1];
	z = normals[i*3+2];
	normtmp = CoordType(x,y,z);
	//normtmp = normtmp/sqrt(normtmp.dot(normtmp));
	(*vi).N() = normtmp;
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
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);//very important !!!
    tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(refmesh);
    float maxDist = m.bbox.Diag();
    float minDist = 1e-10;
    
    vcg::tri::FaceTmark<MyMesh> mf; 
    mf.SetMesh( &m );
    vcg::RayTriangleIntersectionFunctor<true> FintFunct;
    vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
    TriMeshGrid static_grid;    
    static_grid.Set(m.face.begin(), m.face.end());
    // run search 
    for(i=0; i < refmesh.vn; i++)
      {
	hitbool[i] = 0;
	vcg::Ray3f ray;
	Point3f orig = refmesh.vert[i].P();
	Point3f dir = refmesh.vert[i].N();
	Point3f dirOrig = dir;
	
	ray.SetOrigin(orig);
	ray.SetDirection(dir);
	MyFace* f_ptr=GridDoRay(static_grid,FintFunct, mf, ray, maxDist, t);
	
	if (f_ptr)
	  {
	    if (t > 0)
	      {
		MyMesh::CoordType clost = refmesh.vert[i].P()+dir*t;//the hit point
		int f_i = vcg::tri::Index(m, f_ptr);
		MyMesh::CoordType ti = (m.face[f_i].V(0)->N()+m.face[f_i].V(1)->N()+m.face[f_i].V(2)->N())/3;//the smoothed normal at that point
		
		ioclost[i*3] = clost[0];
		ioclost[i*3+1] =clost[1];
		ioclost[i*3+2] =clost[2];
		dis[i]=t;
		
		hitbool[i] = 1;
	      }
	  }
      }
  }
}
	 
	
