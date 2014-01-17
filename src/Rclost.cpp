
#include "typedef.h"
//#include <wrap/ply/plylib.cpp>
#include <cmath>
//#include <utility>  
//#include <algorithm>
extern "C" {

  void Rclost(double *vb ,int *dim, int *it, int *dimit, double *ioclost, int *clostDim, double *normals, double *dis,int *sign,int *border, int *barycentric, double *barycoord,int *faceptr, int* smooth)
  {
    /*typedef MyMesh::CoordType CoordType;
      typedef  MyMesh::ScalarType ScalarType;
    */
    typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid; 
    //typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    ScalarType x,y,z;
    int i;
    MyMesh m;
    MyMesh refmesh;
    MyMesh outmesh;
    MyMesh::CoordType baryco;
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
    SimpleTempData<MyMesh::FaceContainer,int> indices(m.face);
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
    
    VertexIterator vi=m.vert.begin();
    for (i=0; i < d; i++) {
      ivp[i]=&*vi;
      x = vb[i*3];
      y = vb[i*3+1];
      z=  vb[i*3+2];
      (*vi).P() = CoordType(x,y,z);
      ++vi;
    }
    int itx,ity,itz;
    FaceIterator fi=m.face.begin();
    for (i=0; i < faced ; i++) {
      indices[fi] = i;
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
    
    for (i=0; i < dref; i++) {
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
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(m);//very important !!!
    //tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(m);
    if (*smooth == 1)
      tri::Smooth<MyMesh>::VertexNormalLaplacian(m,2,false);
    tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
    float maxDist = m.bbox.Diag()*2;
    float minDist = 1e-10;
    vcg::tri::FaceTmark<MyMesh> mf; 
    mf.SetMesh( &m );
    vcg::face::PointDistanceBaseFunctor<float> PDistFunct;
    TriMeshGrid static_grid;    
    static_grid.Set(m.face.begin(), m.face.end());
    tri::UpdateFlags<MyMesh>::FaceBorderFromNone(m);
    tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(m);
    MyMesh::CoordType tt;
    for(i=0; i < refmesh.vn; i++) {
      border[i]=0;
      Point3f& currp = refmesh.vert[i].P();
      Point3f& clost = outmesh.vert[i].P();
      MyFace* f_ptr= GridClosest(static_grid, PDistFunct, mf, currp, maxDist, minDist, clost);
      if (f_ptr) {
	if ((*f_ptr).IsS())
	  border[i]=1;
	
	faceptr[i] = indices[f_ptr];
	int f_i = vcg::tri::Index(m, f_ptr);
	tt = currp*0;
	
	for (int j=0; j <3;j++) {
	  if (&(m.face[f_i].V(j)->N())) {
	    Point3f vdist = m.face[f_i].V(j)->P() - clost;
	    float weight = sqrt(vdist.dot(vdist));
	    if (weight > 0)
	      weight = 1/weight;
	    else 
	      weight = 1e12;
	    
	    tt +=(m.face[f_i].V(j)->N()*weight);
	  }
	}
	if (*barycentric == 1) {
	  baryco = currp*0;
	  InterpolationParameters<MyFace,ScalarType>(*f_ptr,f_ptr->N(),clost,baryco);
	}
      }
      float vl = sqrt(tt.dot(tt));
      if (vl > 0) {//check for zero length normals
	tt=tt/vl;
      }   
      
      dis[i] = minDist;
      if (signo == 1) {
	Point3f dif = clost - currp;
	float sign = dif.dot(tt);	
	if (sign < 0)
	  dis[i] = -dis[i] ;
      }
      //write back output
      ioclost[i*3] = clost[0];
      ioclost[i*3+1] =clost[1];
      ioclost[i*3+2] =clost[2];
      normals[i*3] = tt[0];
      normals[i*3+1] = tt[1];    
      normals[i*3+2] = tt[2];
      if(*barycentric == 1) {
	barycoord[i*3] = baryco[0];
	barycoord[i*3+1] = baryco[1];
	barycoord[i*3+2] = baryco[2];
      }
    }
  }
}
