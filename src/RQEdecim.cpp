// Author: Stefan Schlager
// This is basically an adaption 
// of tridecimator included in the vcglib sources
// to work with R

#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
// stuff to define the mesh
#include <vcg/complex/complex.h>
#include <vcg/math/quadric.h>
#include <vcg/complex/algorithms/clean.h>

// update
#include <vcg/complex/algorithms/update/topology.h>
//using namespace std;
#include <vcg/complex/algorithms/local_optimization.h>
#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
#include <vcg/container/simple_temporary_data.h>
//#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>

#include "RvcgIO.h"
#include <Rcpp.h>
using namespace vcg;
using namespace tri;
using namespace Rcpp;

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
				vertex::BitFlags  >
{
		 public:
		   vcg::math::Quadric<double> &Qd() {return q;}
		 private:
		   math::Quadric<double> q;
};

class CEdge : public Edge< CUsedTypes> {};

typedef BasicVertexPair<CVertex> VertexPair;

class CFace    : public Face< CUsedTypes,
			      face::VFAdj,
			      //face::FFAdj,
			      face::VertexRef,
			      face::BitFlags > {};

// the main mesh class
class CMeshDec: public vcg::tri::TriMesh<std::vector<CVertex>, 
					 std::vector<CFace> > {};

class CTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< CMeshDec, 
								 VertexPair, 
								 CTriEdgeCollapse,
								 QInfoStandard<CVertex> > {
  public:
  typedef  vcg::tri::TriEdgeCollapseQuadric< CMeshDec,  
					     VertexPair, 
					     CTriEdgeCollapse, 
					     QInfoStandard<CVertex>  > TECQ;
  typedef  CMeshDec::VertexType::EdgeType EdgeType;
  inline CTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
  
};
typedef CMeshDec::VertexIterator VertexIterator;
typedef CMeshDec::FacePointer  FacePointer;
typedef CMeshDec::FaceIterator   FaceIterator;
typedef CMeshDec::CoordType CoordType;
typedef CMeshDec::ScalarType ScalarType;
typedef CMeshDec::VertexPointer VertexPointer;



RcppExport SEXP RQEdecim(SEXP _vb , SEXP _it, SEXP _Finsize, SEXP _boolparams, SEXP _doubleparams)
{
  // declare Mesh and helper variables
  int i;
  CMeshDec m;
  VertexIterator vi;
  FaceIterator fi;
    
  int check = Rvcg::IOMesh<CMeshDec>::RvcgReadR(m,_vb,_it);
  if (check == 1) {
    Rprintf("%s\n","Warning: mesh has no faces, nothing done");
    return Rcpp::List::create(Rcpp::Named("vb") = _vb,
			    Rcpp::Named("normals") = 0,
			    Rcpp::Named("it") = _it
			    );
  }  else {
  Rcpp::LogicalVector boolparams(_boolparams); 
  Rcpp::NumericVector doubleparams(_doubleparams);
  //boolparams: 0=topo 1=quality 2=boundary 3=optiplace 4=scaleindi, 5=  normcheck, 6=safeheap)
  //doubleparams: 0 = qthresh, 1 = boundweight 2=normalthr
  int FinalSize = Rcpp::as<int>(_Finsize);
  
  //initiate decimation process
  TriEdgeCollapseQuadricParameter qparams;
  float TargetError=std::numeric_limits<float>::max();
  qparams.QualityThr = doubleparams(0);
  qparams.BoundaryWeight = doubleparams(1);
  qparams.NormalThrRad = doubleparams(2);

  qparams.PreserveTopology = boolparams(0);
  qparams.QualityCheck = boolparams(1);
  qparams.PreserveBoundary = boolparams(2);
  qparams.OptimalPlacement = boolparams(3);
  qparams.ScaleIndependent = boolparams(4);
  qparams.NormalCheck = boolparams(5);
  qparams.SafeHeapUpdate = boolparams(6);
   
  tri::Clean<CMeshDec>::RemoveDuplicateVertex(m);
  tri::Clean<CMeshDec>::RemoveUnreferencedVertex(m);
  
  Rprintf("reducing it to %i faces\n",FinalSize);
    
  vcg::tri::UpdateBounding<CMeshDec>::Box(m);
    
  // decimator initialization
  vcg::LocalOptimization<CMeshDec> DeciSession(m,&qparams);
    
  DeciSession.Init<CTriEdgeCollapse>();
  Rprintf("Initial Heap Size %i\n",int(DeciSession.h.size()));
    
  DeciSession.SetTargetSimplices(FinalSize);
  DeciSession.SetTimeBudget(0.5f);
  if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);
    
  while(DeciSession.DoOptimization() && m.fn>FinalSize && DeciSession.currMetric < TargetError){}
    
  vcg::tri::Allocator< CMeshDec >::CompactVertexVector(m);
  vcg::tri::Allocator< CMeshDec >::CompactFaceVector(m);
  SimpleTempData<CMeshDec::VertContainer,int> indices(m.vert);
  tri::UpdateNormal<CMeshDec>::PerVertexAngleWeighted(m);
  tri::UpdateNormal<CMeshDec>::NormalizePerVertex(m);
  Rprintf("Result: %d vertices and %d faces.\nEstimated error: %g \n",m.vn,m.fn,DeciSession.currMetric);
    
  Rcpp::NumericMatrix vb(3, m.vn), normals(3, m.vn);
  Rcpp::IntegerMatrix itout(3, m.fn);
  
  //write back data
  vi=m.vert.begin();
  for (i=0;  i < m.vn; i++) {
    indices[vi] = i;//important: updates vertex indices
    vb(0,i) = (*vi).P()[0];
    vb(1,i) = (*vi).P()[1];
    vb(2,i) = (*vi).P()[2];
    normals(0,i) = (*vi).N()[0];
    normals(1,i) = (*vi).N()[1];
    normals(2,i) = (*vi).N()[2];
    ++vi;
  }
  
  FacePointer fp;
  fi=m.face.begin();
  for (i=0; i < m.fn;i++) {
    fp=&(*fi);
    if( ! fp->IsD() ) {
      itout(0,i) = indices[fp->cV(0)]+1;
      itout(1,i) = indices[fp->cV(1)]+1;
      itout(2,i) = indices[fp->cV(2)]+1;
      ++fi;
    }
  }
  return Rcpp::List::create(Rcpp::Named("vb") = vb,
			    Rcpp::Named("normals") = normals,
			    Rcpp::Named("it") = itout
			    );
    }
}
   

 
  
