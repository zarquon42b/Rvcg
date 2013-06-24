#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
// stuff to define the mesh
#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/math/quadric.h>
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

RcppExport  SEXP RgetEdge(SEXP _vb, SEXP _it, SEXP _unique)
{
  int i, j;
  CMesh1 m;
  VertexIterator vi;
  FaceIterator fi;
  bool unique = Rcpp::as<bool>(_unique);  
 
  // allocate and fill mesh
  Rvcg::IOMesh<CMesh1>::RvcgReadR(m,_vb,_it);
  
  // create int indices per face and per vertex to return to R
  SimpleTempData<CMesh1::VertContainer,int> indices(m.vert);
  SimpleTempData<CMesh1::FaceContainer,int> indicesf(m.face);
  vi = m.vert.begin();
  for (i=0; i < m.vn; i++) 
    {indices[vi] = i;
      ++vi;
    }
  fi = m.face.begin();
  for (i=0; i < m.fn ; i++) 
    {indicesf[fi] = i;
      ++fi;
    }
    
  std::vector<SimpleEdge> Edges;
  typename std::vector< SimpleEdge >::iterator ei;
  typename std::vector< SimpleEdge >::size_type size;
  tri::UpdateFlags<CMesh1>::VertexBorderFromNone(m);
  tri::UpdateSelection<CMesh1>::VertexFromBorderFlag(m);
  tri::UpdateTopology<CMesh1>::FaceFace(m);
  tri::UpdateFlags<CMesh1>::FaceBorderFromNone(m); 
    
  if (unique)
    tri::UpdateTopology<CMesh1>::FillUniqueEdgeVector(m,Edges,true);
  else
    tri::UpdateTopology<CMesh1>::FillEdgeVector(m,Edges,true);
   
  size=Edges.size();
  Rcpp::IntegerVector facept(size), border(size);
  Rcpp::IntegerMatrix edges(size,2);
  border = border * 0;
  EdgePointer ep;
  VertexPointer vp , vp1;
  FacePointer fp;
  i=0;
    
  // for(ei=Edges.begin(); ei!=Edges.end(); ++ei)
  for (i = 0;i < size;i++)
    {
      vp=Edges[i].v[0];
      vp1=Edges[i].v[1];
      fp=Edges[i].f;
      if( (*fp).IsB(Edges[i].z)==true )
	border[i] = 1;
		
      edges(i,0)=indices[vp]+1;
      edges(i,1)=indices[vp1]+1;
      facept[i] = indicesf[Edges[i].f]+1;
    }
    
  return Rcpp::List::create(Rcpp::Named("edges") = edges,
			    Rcpp::Named("facept") = facept,    
			    Rcpp::Named("border") = border
			    );
}
    
