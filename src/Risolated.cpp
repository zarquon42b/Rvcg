// Author: Stefan Schlager
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

using namespace vcg;
using namespace tri;
// The class prototypes.
class CVertex;
class CEdge;
class CFace;

struct CUsedTypes: public UsedTypes<Use<CVertex>::AsVertexType,Use<CEdge>::AsEdgeType,Use<CFace>::AsFaceType>{};

class CVertex  : public Vertex< CUsedTypes,
  vertex::VFAdj,
  vertex::Coord3f,
  vertex::Normal3f,
  vertex::Mark,
				vertex::BitFlags  >{};

class CEdge : public Edge< CUsedTypes> {};


class CFace    : public Face< CUsedTypes,
  face::VFAdj,
  face::VertexRef,
  face::FFAdj,
  face::Mark,
  face::BitFlags > {};

// the main mesh class
class CMesh1    : public vcg::tri::TriMesh<std::vector<CVertex>, std::vector<CFace> > {};
 typedef CMesh1::VertexIterator VertexIterator;
    typedef CMesh1::FacePointer  FacePointer;
    typedef CMesh1::FaceIterator   FaceIterator;
    typedef CMesh1::CoordType CoordType;
    typedef CMesh1::ScalarType ScalarType;
    typedef CMesh1::VertexPointer VertexPointer;
typedef Point3<CMesh1::ScalarType> Point3x;
typedef std::vector<Point3x> Hole;


typedef CMesh1::FaceContainer FaceContainer;


extern "C" {
  void Risolated(double *vb ,int *dim, int *it, int *dimit,double *diam, double *normals,int *facenum)
  {
    // typedefs
   
    ScalarType x,y,z;
    int i;
    CMesh1 m;
    int d = *dim;
    int faced = *dimit;
    float diameter = *diam;
    int connect = *facenum;
    //create mesh
    vcg::tri::Allocator<CMesh1>::AddVertices(m,d);
    vcg::tri::Allocator<CMesh1>::AddFaces(m,faced);
    typedef CMesh1::VertexPointer VertexPointer;
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
    //tri::Clean<CMesh1>::RemoveDuplicateVertex(m);
    tri::UpdateTopology<CMesh1>::FaceFace(m);
    tri::UpdateTopology<CMesh1>::VertexFace(m);
    std::pair<int,int> delInfo;
      std::vector< std::pair<int,CMesh1::FacePointer> > CCV;
    int TotalCC=tri::Clean<CMesh1>::ConnectedComponents(m, CCV);
    //printf("%i\n",TotalCC);
    std::vector<float> chunks;
    std::vector<int> chunkface;
    //int CCm = tri::Clean<CMesh1>::ConnectedComponents(m);
    tri::ConnectedIterator<CMesh1> ci;
    for(unsigned int i1=0;i1<CCV.size();++i1)
      {
	Box3f bb;
	std::vector<CMesh1::FacePointer> FPV;
	for(ci.start(m,CCV[i1].second);!ci.completed();++ci)
	  {
	    FPV.push_back(*ci);
	    bb.Add((*ci)->P(0));
	    bb.Add((*ci)->P(1));
	    bb.Add((*ci)->P(2));
	  } 
	float diag = bb.Diag();
	chunks.push_back(diag);
	chunkface.push_back(CCV[i1].first);
      }
    // int myints[] = {3,7,2,5,6,4,9};
    //int ppp = *std::max_element(myints,myints+7);
    if (diameter == 0)
      diameter = *std::max_element(chunks.begin(),chunks.end());
    
    if (connect < 0)
      delInfo= tri::Clean<CMesh1>::RemoveSmallConnectedComponentsDiameter(m,diameter);
     
    else 
//if (connect == 0)
      {
	if (connect == 0)
	  connect = *std::max_element(chunkface.begin(),chunkface.end());
	delInfo = tri::Clean<CMesh1>::RemoveSmallConnectedComponentsSize(m,connect);
      }
      
    printf("Removed %i connected components out of %i\n", delInfo.second, delInfo.first); 
	//int unref =  tri::Clean<CMesh1>::RemoveUnreferencedVertex(m);
    
      
	
    vcg::tri::Allocator< CMesh1 >::CompactVertexVector(m);
    vcg::tri::Allocator< CMesh1 >::CompactFaceVector(m);
    SimpleTempData<CMesh1::VertContainer,int> indices(m.vert);
    tri::UpdateNormals<CMesh1>::PerVertexAngleWeighted(m);
    tri::UpdateNormals<CMesh1>::NormalizeVertex(m);

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



