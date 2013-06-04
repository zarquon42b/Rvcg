// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
// stuff to define the mesh
#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/edge/base.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/edges.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/clean.h>
// io
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>
#include <wrap/io_trimesh/export_ply.h>
#include <wrap/io_trimesh/export_ply.h>
//#include <wrap/ply/plylib.cpp>
// update
#include <vcg/complex/algorithms/update/topology.h>
//using namespace std;
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <vcg/container/simple_temporary_data.h>
#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <vcg/simplex/face/pos.h>

#include <iostream>
using namespace vcg;
using namespace tri;
// The class prototypes.
class CVertex1;
class CEdge1;
class CFace1;

struct CUsedTypes1: public UsedTypes<Use<CVertex1>::AsVertexType,Use<CEdge1>::AsEdgeType,Use<CFace1>::AsFaceType>{};

class CVertex1  : public Vertex< CUsedTypes1,
				 vertex::VFAdj,
				 vertex::Coord3f,
				 vertex::Normal3f,
				 vertex::Mark,
				 vertex::BitFlags  >{};

class CEdge1 : public Edge< CUsedTypes1> {};


class CFace1    : public Face< CUsedTypes1,
			       face::VFAdj,
			       face::VertexRef,
			       face::FFAdj,
			       face::Mark,
			       face::BitFlags > {};

// the main mesh class
class CMesh1    : public vcg::tri::TriMesh<std::vector<CVertex1>, std::vector<CFace1> > {};
typedef CMesh1::VertexIterator VertexIterator;
typedef CMesh1::FacePointer  FacePointer;
typedef CMesh1::FaceIterator   FaceIterator;
typedef CMesh1::EdgePointer   EdgePointer;
typedef CMesh1::EdgeIterator   EdgeIterator;
typedef CMesh1::CoordType CoordType;
typedef CMesh1::ScalarType ScalarType;
typedef CMesh1::VertexPointer VertexPointer;
typedef Point3<CMesh1::ScalarType> Point3x;
typedef std::vector<Point3x> Hole;
//typedef face::Pos<CFace1>    PosType;
//typedef vcg::face::Pos<CFace1> PsType;

typedef CMesh1::FaceContainer FaceContainer;

extern "C" {

  void RgetEdge(double *vb ,int *dim, int *it, int *dimit, int *edgecount, int *edges,int *facept, int *border,int *unique)
  {
    /*typedef MyMesh::CoordType CoordType;
      typedef  MyMesh::ScalarType ScalarType;
    
      typedef vcg::SpatialHashTable<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;*/ 
    //typedef vcg::GridStaticPtr<MyMesh::FaceType, MyMesh::ScalarType> TriMeshGrid;
    ScalarType x,y,z;
    int i;
    
    CMesh1 m;
    *border=*border*0;
    // section read from input
    const int d = *dim;
    const int faced = *dimit;
   
   
   
    //Allocate mesh

    typedef UpdateTopology<CMesh1>::PEdge SimpleEdge;
    vcg::tri::Allocator<CMesh1>::AddVertices(m,d);
    vcg::tri::Allocator<CMesh1>::AddFaces(m,faced);
    typedef CMesh1::VertexPointer VertexPointer;
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
    SimpleTempData<CMesh1::VertContainer,int> indices(m.vert);
    SimpleTempData<CMesh1::FaceContainer,int> indicesf(m.face);

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
    std::vector<SimpleEdge> Edges;
    typename std::vector< SimpleEdge >::iterator ei;
    typename std::vector< SimpleEdge >::size_type size;
    tri::UpdateFlags<CMesh1>::VertexBorderFromNone(m);
    tri::UpdateSelection<CMesh1>::VertexFromBorderFlag(m);
    tri::UpdateTopology<CMesh1>::FaceFace(m);
    tri::UpdateFlags<CMesh1>::FaceBorderFromNone(m); 
    if (*unique == 1)
      tri::UpdateTopology<CMesh1>::FillUniqueEdgeVector(m,Edges,true);
    else
      tri::UpdateTopology<CMesh1>::FillEdgeVector(m,Edges,true);
    EdgePointer ep;
    VertexPointer vp , vp1;
    FacePointer fp;
    i=0;
    size=Edges.size();
    // for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
    for (i=0;i<size;i++)
      {
	vp=Edges[i].v[0];
	vp1=Edges[i].v[1];
	//if ((*vp).IsS() &&(*vp1).IsS())
	//  border[i] = 1;
	fp=Edges[i].f;
	if( (*fp).IsB(Edges[i].z)==true )
	  border[i] = 1;
	edges[i*2]=indices[vp];
	edges[i*2+1]=indices[vp1];
	facept[i] = indicesf[Edges[i].f];
	//printf("%d\n",facept[i]);	
      }
   
    *edgecount= size;
  }

}
    
