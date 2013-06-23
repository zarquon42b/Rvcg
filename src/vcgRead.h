#include <vector>
#include <limits>
#include <Rcpp.h>
#include <vcg/complex/complex.h>
#include<vcg/complex/allocate.h>
#include <wrap/callback.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/space/point3.h>

//using namespace vcg;
//using namespace tri;

//using namespace Rcpp;
//using namespace std;

template <class ReadMeshType>
class ReadMesh
{

public:
  typedef ReadMeshType MeshType; 
  typedef typename MeshType::CoordType CoordType;
  typedef typename MeshType::ScalarType     ScalarType;
  typedef typename MeshType::VertexType     VertexType;
  typedef typename MeshType::VertexPointer  VertexPointer;
  typedef typename MeshType::VertexIterator VertexIterator;
  typedef typename MeshType::FaceType       FaceType;
  typedef typename MeshType::FacePointer    FacePointer;
  typedef typename MeshType::FaceIterator   FaceIterator;
  typedef typename MeshType::FaceContainer  FaceContainer;
  typedef typename MeshType::VertContainer  VertContainer;
  

static void vcgReadR(MeshType &m, SEXP _vb, SEXP _it)
{
  Rcpp::IntegerMatrix it(_it);
  Rcpp::NumericMatrix vb(_vb);
  int d =  vb.ncol();
  int faced = it.ncol();
  ScalarType x,y,z;
  int i, j;
  vcg::tri::Allocator<MeshType>::AddVertices(m,d);
  vcg::tri::Allocator<MeshType>::AddFaces(m,faced);
  std::vector<VertexPointer> ivp;
  ivp.resize(d);
  
  vcg::SimpleTempData<typename MeshType::FaceContainer, int> indicesf(m.face);
  vcg::SimpleTempData<typename MeshType::VertContainer, int> indices(m.vert);
  VertexIterator vi = m.vert.begin();
  for (i=0; i < d; i++) 
    {
      ivp[i]=&*vi;
      x = vb(0,i);
      y = vb(1,i);
      z = vb(2,i);
      (*vi).P() = CoordType(x,y,z);
      ++vi;
    } 
    
  int itx,ity,itz;
  FaceIterator fi=m.face.begin();
  for (i=0; i < faced ; i++) 
    {
      indicesf[fi] = i;
      itx = it(0,i);
      ity = it(1,i);
      itz = it(2,i);
      (*fi).V(0)=ivp[itx];
      (*fi).V(1)=ivp[ity];
      (*fi).V(2)=ivp[itz];
      ++fi;
    }
  //return(0);
  
}
  

};
