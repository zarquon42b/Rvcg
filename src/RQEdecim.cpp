// Author: Stefan Schlager
// This is basically an adaption 
// of tridecimator included in the vcglib sources
// to work with R

#include <vector>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
// stuff to define the mesh
/*#include <vcg/simplex/vertex/base.h>
#include <vcg/simplex/face/base.h>
#include <vcg/simplex/edge/base.h>*/
#include <vcg/complex/complex.h>
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
				vertex::BitFlags  >{
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
class CMeshDec    : public vcg::tri::TriMesh<std::vector<CVertex>, std::vector<CFace> > {};

class CTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< CMeshDec, VertexPair, CTriEdgeCollapse, QInfoStandard<CVertex>  > {
public:
  typedef  vcg::tri::TriEdgeCollapseQuadric< CMeshDec,  VertexPair, CTriEdgeCollapse, QInfoStandard<CVertex>  > TECQ;
  typedef  CMeshDec::VertexType::EdgeType EdgeType;
  inline CTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}

};

extern "C" {

  void RQEdecim(double *vb ,int *dim, int *it, int *dimit,int *Finsize,double *normals,int *topo,int *qual,int *bound)
  {
    // typedefs
    typedef CMeshDec::VertexIterator VertexIterator;
    typedef CMeshDec::FacePointer  FacePointer;
    typedef CMeshDec::FaceIterator   FaceIterator;
    typedef CMeshDec::CoordType CoordType;
    typedef CMeshDec::ScalarType ScalarType;
    typedef CMeshDec::VertexPointer VertexPointer;
    //set local variables
    ScalarType x,y,z;
    int i;
    int FinalSize=*Finsize;
    int d = *dim;
    int faced = *dimit;
    CMeshDec m;
   
    // fill mesh with data from R workspace
    vcg::tri::Allocator<CMeshDec>::AddVertices(m,d);
    vcg::tri::Allocator<CMeshDec>::AddFaces(m,faced);
    //VertexPointer ivp[d];
    std::vector<VertexPointer> ivp;
    ivp.resize(d);
    VertexIterator vi=m.vert.begin();
    for (i=0; i<d; i++) 
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
    for (i=0; i < faced; i++) 
      {
	itx = it[i*3];
	ity = it[i*3+1];
	itz = it[i*3+2];
	(*fi).V(0)=ivp[itx];
	(*fi).V(1)=ivp[ity];
	(*fi).V(2)=ivp[itz];
	++fi;
      }
    //tri::io::ExporterPLY<CMeshDec>::Save(m,"tt.ply",tri::io::Mask::IOM_VERTNORMAL, false); // in ASCII
    //initiate decimation process
    TriEdgeCollapseQuadricParameter qparams;
    float TargetError=std::numeric_limits<float>::max();
    ///qparams.QualityThr =.3;
    qparams.QualityCheck = true;
    qparams.OptimalPlacement = true;
    qparams.PreserveTopology	= true;
    qparams.PreserveBoundary	= true;
    if (*topo == 0)
      qparams.PreserveTopology	= false;
    if (*qual == 0)
      qparams.OptimalPlacement = true;
    if (*bound == 0)
      qparams.PreserveBoundary	= false;
   
   
   
    int dup = tri::Clean<CMeshDec>::RemoveDuplicateVertex(m);
    int unref =  tri::Clean<CMeshDec>::RemoveUnreferencedVertex(m);
   
    
    printf("reducing it to %i faces\n",FinalSize);
    
    vcg::tri::UpdateBounding<CMeshDec>::Box(m);
    
    // decimator initialization
    vcg::LocalOptimization<CMeshDec> DeciSession(m,&qparams);
    
    int t1=clock();
    DeciSession.Init<CTriEdgeCollapse>();
    int t2=clock();
    //printf("Initial Heap Size %i\n",int(DeciSession.h.size()));
    
    DeciSession.SetTargetSimplices(FinalSize);
    DeciSession.SetTimeBudget(0.5f);
    if(TargetError< std::numeric_limits<float>::max() ) DeciSession.SetTargetMetric(TargetError);
    
    while(DeciSession.DoOptimization() && m.fn>FinalSize && DeciSession.currMetric < TargetError)
      // printf("Final Mesh size: %7i heap sz %9i err %9g \r",m.fn, int(DeciSession.h.size()),DeciSession.currMetric);
    
      int t3=clock();
   
    //printf("\nCompleted in (%i+%i) msec\n",t2-t1,t3-t2);
    
    vcg::tri::Allocator< CMeshDec >::CompactVertexVector(m);
    vcg::tri::Allocator< CMeshDec >::CompactFaceVector(m);
    SimpleTempData<CMeshDec::VertContainer,int> indices(m.vert);
    tri::UpdateNormal<CMeshDec>::PerVertexAngleWeighted(m);
    tri::UpdateNormal<CMeshDec>::NormalizePerVertex(m);
    printf("Result: %d vertices and %d faces.\nEstimated error: %g \n",m.vn,m.fn,DeciSession.currMetric);
    
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
    //  printf("%i %i\n",m.vn,m.fn);
    
    
  }
   
}
 
  
