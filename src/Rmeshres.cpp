// stuff to define the mesh
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
//#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/clean.h>
// update
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/container/simple_temporary_data.h>
#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <vcg/simplex/face/pos.h>

#include <../RvcgIO.h>
#include <Rcpp.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;

// The class prototypes.
class CVertex1;
class CEdge1;
class CFace1;

struct CUsedTypes1: public UsedTypes<Use<CVertex1>::AsVertexType,
				     Use<CEdge1>::AsEdgeType,
				     Use<CFace1>::AsFaceType>{};

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
class CMesh1    : public vcg::tri::TriMesh<std::vector<CVertex1>, 
					   std::vector<CFace1> > {};
typedef CMesh1::VertexIterator VertexIterator;
typedef CMesh1::FacePointer  FacePointer;
typedef CMesh1::FaceIterator   FaceIterator;
typedef CMesh1::EdgePointer   EdgePointer;
typedef CMesh1::EdgeIterator   EdgeIterator;
typedef CMesh1::CoordType CoordType;
typedef CMesh1::ScalarType ScalarType;
typedef CMesh1::VertexPointer VertexPointer;
typedef Point3<CMesh1::ScalarType> Point3x;
//typedef std::vector<Point3x> Hole;
typedef CMesh1::FaceContainer FaceContainer;
typedef UpdateTopology<CMesh1>::PEdge SimpleEdge;
  
RcppExport SEXP Rmeshres(SEXP _vb , SEXP _it)
  {
    // declare Mesh and helper variables
    int i, j;
    CMesh1 m;
    VertexIterator vi;
    FaceIterator fi;
    
    Rvcg::IOMesh<CMesh1>::RvcgReadR(m,_vb,_it);
    std::vector<SimpleEdge> Edges;
    typename std::vector< SimpleEdge >::iterator ei;
    typename std::vector< SimpleEdge >::size_type size;
    tri::UpdateTopology<CMesh1>::FaceFace(m);
    tri::UpdateTopology<CMesh1>::FillUniqueEdgeVector(m,Edges,true);
    size=Edges.size();
    double res = 0;
    Point3f tmp0;
    VertexPointer vp , vp1;
    for (i = 0;i < size;i++)
    {
      vp=Edges[i].v[0];
      vp1=Edges[i].v[1];
      tmp0 = vp->P()-vp1->P();
      res = res + sqrt(tmp0.dot(tmp0));
    }
    res = res/size;
    return(wrap(res));
  }

