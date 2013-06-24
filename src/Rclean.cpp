// This is based on code from 
// of trimeshinfo
// included in the vcglib sources
// to work with R
#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
// stuff to define the mesh

#include <vcg/complex/complex.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/flag.h>
#include <vcg/complex/algorithms/clean.h>
// update
#include <vcg/complex/algorithms/update/topology.h>
//using namespace std;

#include <vcg/container/simple_temporary_data.h>
#include<vcg/complex/allocate.h>
//#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <Rcpp.h>
#include <../RvcgIO.h>

using namespace vcg;
using namespace tri;
using namespace Rcpp;
using namespace std;

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
                               face::Normal3f,
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

RcppExport SEXP Rclean(SEXP _vb, SEXP _it, SEXP _type)
{
  // declare Mesh and helper variables
  int select = Rcpp::as<int>(_type);  
  int i, j;
  CMesh2 m;
  VertexIterator vi;
  FaceIterator fi;
  // allocate mesh and fill it
  Rvcg::IOMesh<CMesh2>::RvcgReadR(m,_vb,_it);
  
  // General cleaning and update of topology
  //tri::UpdateFlags<CMesh2>::VertexBorderFromNone(m);
  //tri::UpdateSelection<CMesh2>::VertexFromBorderFlag(m);
  int dupvb = tri::Clean<CMesh2>::RemoveDuplicateVertex(m);
  int dupit = tri::Clean<CMesh2>::RemoveDuplicateFace(m);
  tri::UpdateTopology<CMesh2>::FaceFace(m);
  tri::UpdateTopology<CMesh2>::VertexFace(m);
  vcg::tri::UpdateFlags<CMesh2>::FaceBorderFromFF(m);
  vcg::tri::UpdateFlags<CMesh2>::VertexBorderFromFace(m);
  printf("removed %d duplicate faces and %d duplicate vertices\n",dupit,dupvb);
  //tri::UpdateFlags<CMesh2>::FaceBorderFromNone(m); 
   
  // do all the cleaning
    
  if (select == 1)
    { int unref = tri::Clean<CMesh2>::RemoveUnreferencedVertex(m);
      printf("removed %d unreferenced vertices\n",unref);
    }
  int rem;
  if (select == 2)
    { rem = tri::Clean<CMesh2>::RemoveNonManifoldFace(m);
      printf("removed %d Non-manifold faces\n",rem);
    }
  if (select == 3)
    { rem = tri::Clean<CMesh2>::RemoveDegenerateFace(m);
      printf("removed %d degenerate faces\n",rem);
    }
  if (select == 4)
    { rem = tri::Clean<CMesh2>::RemoveNonManifoldVertex(m);
      printf("removed %d Non-manifold vertices\n",rem);
    }
  if (select == 5)
    { int split =tri::Clean<CMesh2>::SplitNonManifoldVertex(m,0.1);
      printf("split %d non-manifold vertices\n",split);
    }
  
 
  // write back
  vcg::tri::Allocator< CMesh2 >::CompactVertexVector(m);
  vcg::tri::Allocator< CMesh2 >::CompactFaceVector(m);
  tri::UpdateNormal<CMesh2>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<CMesh2>::NormalizePerVertex(m);
  Rcpp::NumericMatrix vbout(3,m.vn), normals(3,m.vn);
  Rcpp::IntegerMatrix itout(3,m.fn);
  
  //write back
 
  vi=m.vert.begin();
  SimpleTempData<CMesh2::VertContainer,int> indiceout(m.vert);
 
  for (i=0;  i < m.vn; i++) 
    {
      if( ! vi->IsD() )
	{
	  indiceout[vi] = i;//important: updates vertex indices
	  vbout(0,i) = (*vi).P()[0];
	  vbout(1,i) = (*vi).P()[1];
	  vbout(2,i) = (*vi).P()[2];
	  normals(0,i) = (*vi).N()[0];
	  normals(1,i) = (*vi).N()[1];
	  normals(2,i) = (*vi).N()[2];
	}
	  ++vi;
    }
  
  FacePointer fp;
  fi=m.face.begin();
       
  for (i=0; i < m.fn; i++) 
    {
      fp=&(*fi);
      if( ! fp->IsD() )
	{
	  itout(0,i) = indiceout[fp->cV(0)]+1;
	  itout(1,i) = indiceout[fp->cV(1)]+1;
	  itout(2,i) = indiceout[fp->cV(2)]+1;
	  ++fi;
	}
    }
  
  return Rcpp::List::create(Rcpp::Named("vb") = vbout,
			    Rcpp::Named("it") = itout,
			    Rcpp::Named("normals") = normals
			    );
}
 

    
