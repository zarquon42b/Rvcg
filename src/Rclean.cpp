// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
// stuff to define the mesh
//#include <vcg/simplex/vertex/base.h>
//#include <vcg/simplex/face/base.h>
//#include <vcg/simplex/edge/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
//#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
//#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/clean.h>
// update
#include <vcg/complex/algorithms/update/topology.h>
//using namespace std;
//#include <vcg/complex/algorithms/local_optimization.h>
//#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <vcg/container/simple_temporary_data.h>
#include<vcg/complex/allocate.h>
//#include <wrap/callback.h>
#include <vcg/complex/append.h>
//#include <vcg/simplex/face/pos.h>
//#include <iostream>
using namespace vcg;
using namespace tri;

// The class prototypes.
class CVertex2;
class CEdge2;
class CFace2;

struct CUsedTypes2: public UsedTypes<Use<CVertex2>::AsVertexType,Use<CEdge2>::AsEdgeType,Use<CFace2>::AsFaceType>{};

class CVertex2  : public Vertex< CUsedTypes2,
				 vertex::VFAdj,
				 vertex::Coord3f,
				 vertex::Normal3f,
				 vertex::Mark,
				 vertex::BitFlags  >{};

class CEdge2 : public Edge< CUsedTypes2> {};


class CFace2    : public Face< CUsedTypes2,
			       face::VFAdj,
			       face::VertexRef,
			       face::FFAdj,
			       face::Mark,
			       face::BitFlags > {};

// the main mesh class
class CMesh2    : public vcg::tri::TriMesh<std::vector<CVertex2>, std::vector<CFace2> > {};
typedef CMesh2::VertexIterator VertexIterator;
typedef CMesh2::FacePointer  FacePointer;
typedef CMesh2::FaceIterator   FaceIterator;
typedef CMesh2::EdgePointer   EdgePointer;
typedef CMesh2::EdgeIterator   EdgeIterator;
typedef CMesh2::CoordType CoordType;
typedef CMesh2::ScalarType ScalarType;
typedef CMesh2::VertexPointer VertexPointer;
typedef Point3<CMesh2::ScalarType> Point3x;
typedef std::vector<Point3x> Hole;


typedef CMesh2::FaceContainer FaceContainer;

extern "C" {

  void Rclean(double *vb ,int *dim, int *it, int *dimit,double *normals)
  {
    ScalarType x,y,z;
    int i;
    
    CMesh2 m;
    
    // section read from input
    int d = *dim;
    int faced = *dimit;
   
   
   
    //Allocate mesh

    typedef UpdateTopology<CMesh2>::PEdge SimpleEdge;
    vcg::tri::Allocator<CMesh2>::AddVertices(m,d);
    vcg::tri::Allocator<CMesh2>::AddFaces(m,faced);
    typedef CMesh2::VertexPointer VertexPointer;
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
    SimpleTempData<CMesh2::VertContainer,int> indices(m.vert);
    SimpleTempData<CMesh2::FaceContainer,int> indicesf(m.face);

    VertexIterator vi=m.vert.begin();
    for (i=0; i < d; i++) 
      {
	indices[vi] = i;
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
	indicesf[fi] = i;
	itx = it[i*3];
	ity = it[i*3+1];
	itz = it[i*3+2];
	(*fi).V(0)=ivp[itx];
	(*fi).V(1)=ivp[ity];
	(*fi).V(2)=ivp[itz];
	++fi;
      }
    
    tri::UpdateFlags<CMesh2>::VertexBorderFromNone(m);
    tri::UpdateSelection<CMesh2>::VertexFromBorderFlag(m);
    tri::UpdateTopology<CMesh2>::FaceFace(m);
    tri::UpdateTopology<CMesh2>::VertexFace(m);
    //tri::UpdateFlags<CMesh2>::FaceBorderFromNone(m); 
   
    // do all the cleaning
    //tri::Clean<CMesh2>::SplitNonManifoldVertex(m,0.1);
    tri::Clean<CMesh2>::RemoveNonManifoldFace(m);
    tri::Clean<CMesh2>::RemoveDegenerateFace(m);
    tri::Clean<CMesh2>::RemoveDuplicateVertex(m);
    tri::Clean<CMesh2>::RemoveDuplicateFace(m);
    tri::Clean<CMesh2>::RemoveUnreferencedVertex(m);
    tri::Clean<CMesh2>::RemoveNonManifoldVertex(m);
    

   

    vcg::tri::Allocator< CMesh2 >::CompactVertexVector(m);
    vcg::tri::Allocator< CMesh2 >::CompactFaceVector(m);
    //SimpleTempData<CMesh2::VertContainer,int> indices(m.vert);
    tri::UpdateNormal<CMesh2>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<CMesh2>::NormalizePerVertex(m);
    i=0;
    //size=Edges.size();
    // for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    
    vi=m.vert.begin();
    for (i=0;  i < m.vn; i++) 
      {
	indices[vi] = i;//important: updates vertex indices
	vb[i*3] = (*vi).P()[0];
	vb[i*3+1] = (*vi).P()[1];
	vb[i*3+2] = (*vi).P()[2];
	normals[i*3] = (*vi).N()[0];
	normals[i*3+1] = (*vi).N()[1];
	normals[i*3+2] = (*vi).N()[2];
	++vi;
      }
    
    FacePointer fp;
    int vv[3];
    *dim = m.vn;
    fi=m.face.begin();
    faced=m.fn;
    
    for (i=0; i < faced;i++) 
      {
	fp=&(*fi);
	if( ! fp->IsD() )
	  {
	    vv[0]=indices[fp->cV(0)];
	    vv[1]=indices[fp->cV(1)];
	    vv[2]=indices[fp->cV(2)];
	    it[i*3]=vv[0];
	    it[i*3+1]=vv[1];
	    it[i*3+2]=vv[2];
	    ++fi;
	  }
      }
    *dimit=m.fn;
  }

}
    
