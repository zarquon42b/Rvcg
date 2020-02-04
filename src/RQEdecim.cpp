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
#include <RcppArmadillo.h>
#ifndef RcppExport
#define RcppExport extern "C"
#endif


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



RcppExport SEXP RQEdecim(SEXP mesh_, SEXP Finsize_, SEXP boolparams_, SEXP doubleparams_,SEXP silent_)
{
  try {
    // declare Mesh and helper variables
    int i;
    CMeshDec m;
    VertexIterator vi;
    FaceIterator fi;
    bool silent = as<bool>(silent_);
    int check = Rvcg::IOMesh<CMeshDec>::mesh3d2Rvcg(m,mesh_);
    if (check == 1) {
      ::Rf_error("mesh has no faces");
    }  else {
      Rcpp::LogicalVector boolparams(boolparams_); 
      Rcpp::NumericVector doubleparams(doubleparams_);
      //boolparams: 0=topo 1=quality 2=boundary 3=optiplace 4=scaleindi, 5=  normcheck, 6=safeheap)
      //doubleparams: 0 = qthresh, 1 = boundweight 2=normalthr
      int FinalSize = Rcpp::as<int>(Finsize_);
  
      //initiate decimation process
      TriEdgeCollapseQuadricParameter qparams;
      float TargetError=std::numeric_limits<float>::max();
      qparams.QualityThr = doubleparams[0];
      qparams.BoundaryQuadricWeight = doubleparams[1];
      qparams.NormalThrRad = doubleparams[2];

      qparams.PreserveTopology = boolparams[0];
      qparams.QualityCheck = boolparams[1];
      qparams.PreserveBoundary = boolparams[2];
      qparams.OptimalPlacement = boolparams[3];
      qparams.ScaleIndependent = boolparams[4];
      qparams.NormalCheck = boolparams[5];
      qparams.QualityWeightFactor = boolparams[6];
   
      tri::Clean<CMeshDec>::RemoveDuplicateVertex(m);
      tri::Clean<CMeshDec>::RemoveUnreferencedVertex(m);
      if (!silent)
	Rprintf("reducing it to %i faces\n",FinalSize);
    
      vcg::tri::UpdateBounding<CMeshDec>::Box(m);
    
      // decimator initialization
      vcg::LocalOptimization<CMeshDec> DeciSession(m,&qparams);
    
      DeciSession.Init<CTriEdgeCollapse>();
      if (!silent)
	Rprintf("Initial Heap Size %i\n",int(DeciSession.h.size()));
    
      DeciSession.SetTargetSimplices(FinalSize);
      DeciSession.SetTimeBudget(0.5f);
      if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);
    
      while(m.fn > FinalSize && DeciSession.currMetric < TargetError){
	DeciSession.DoOptimization();
      }
    
      vcg::tri::Allocator< CMeshDec >::CompactVertexVector(m);
      vcg::tri::Allocator< CMeshDec >::CompactFaceVector(m);
      SimpleTempData<CMeshDec::VertContainer,int> indices(m.vert);
      tri::UpdateNormal<CMeshDec>::PerVertexAngleWeighted(m);
      tri::UpdateNormal<CMeshDec>::NormalizePerVertex(m);
      if (!silent)
	Rprintf("Result: %d vertices and %d faces.\nEstimated error: %g \n",m.vn,m.fn,DeciSession.currMetric);
    
      
      //write back data
      List out = Rvcg::IOMesh<CMeshDec>::RvcgToR(m);
      return out;
    }
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
   

 
  
