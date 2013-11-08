#include "typedef.h"
  
  
extern "C" {

  void Rborder(double *vb ,int *dim, int *it, int *dimit, int *bordervb, int *borderit)
  {
    ScalarType x,y,z;
    int i;
    MyMesh m;
    MyMesh refmesh;
    MyMesh outmesh;
    // section read from input
    const int d = *dim;
    const int faced = *dimit;
    
    vcg::tri::Allocator<MyMesh>::AddVertices(m,d);
    vcg::tri::Allocator<MyMesh>::AddFaces(m,faced);
    typedef MyMesh::VertexPointer VertexPointer;
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
    
    //write back border vertices
    vi=m.vert.begin();
    for(i=0; i < m.vn; i++) {
      bordervb[i]=0;
      if ((*vi).IsS())
	bordervb[i]=1;
      ++vi;    
    }
    fi=m.face.begin();
    for(i=0; i < m.fn; i++) { 
      borderit[i]=0;
      if ((*fi).IsS())
	borderit[i]=1;
      ++fi;    
    }
  }
}
