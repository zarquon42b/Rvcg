#include <vector>
#include <limits>
#include <Rcpp.h>
#include <vcg/complex/complex.h>
//#include <vcg/complex/allocate.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/space/point3.h>

namespace Rvcg
{
  template <class IOMeshType>
    class IOMesh
    {
    public:
      typedef IOMeshType MeshType; 
      typedef typename MeshType::CoordType      CoordType;
      typedef typename MeshType::ScalarType     ScalarType;
      typedef typename MeshType::VertexType     VertexType;
      typedef typename MeshType::VertexPointer  VertexPointer;
      typedef typename MeshType::VertexIterator VertexIterator;
      typedef typename MeshType::FaceType       FaceType;
      typedef typename MeshType::FacePointer    FacePointer;
      typedef typename MeshType::FaceIterator   FaceIterator;
      typedef typename MeshType::FaceContainer  FaceContainer;
      typedef typename MeshType::VertContainer  VertContainer;
  
      // Fill an empty mesh with vertices and faces from R
      static int RvcgReadR(MeshType &m, SEXP vb_, SEXP it_= Rcpp::wrap(0), SEXP normals_ = Rcpp::wrap(0))
      {
	//insert vertices
	if (Rf_isMatrix(vb_) && VertexType::HasCoord() ) {
	  Rcpp::NumericMatrix vb(vb_);
	  int d =  vb.ncol();
	  ScalarType x,y,z;	 
	  vcg::tri::Allocator<MeshType>::AddVertices(m,d);
	  std::vector<VertexPointer> ivp;
	  ivp.resize(d);
	  vcg::SimpleTempData<typename MeshType::VertContainer, int> indices(m.vert);
	  VertexIterator vi = m.vert.begin();
	  for (int i=0; i < d; i++) {
	    ivp[i]=&*vi;
	    x = vb(0,i);
	    y = vb(1,i);
	    z = vb(2,i);
	    (*vi).P() = CoordType(x,y,z);
	    ++vi;
	  }
	  
	  //insert vertex normals
	  if (Rf_isMatrix(normals_) && vcg::tri::HasPerVertexNormal(m)) {
	    Rcpp::NumericMatrix normals(normals_);
	    if (normals.ncol() != d) {
	      Rprintf("number of normals is not equal to number of vertices");
	    } else {
	      vcg::SimpleTempData<typename MeshType::VertContainer, int> indices(m.vert);
	      vi = m.vert.begin();
	      for (int i=0; i < d; i++) {
		ivp[i]=&*vi;
		x = normals(0,i);
		y = normals(1,i);
		z = normals(2,i);
		(*vi).N() = CoordType(x,y,z);
		++vi;
	      }
	    }
	  
	  }
	  
	  //process faces but check attributes and input first
	  if (Rf_isMatrix(it_) && FaceType::HasVertexRef()) {
	    Rcpp::IntegerMatrix it(it_);
	    int faced = it.ncol();
	    vcg::tri::Allocator<MeshType>::AddFaces(m,faced);
	    vcg::SimpleTempData<typename MeshType::FaceContainer, int> indicesf(m.face);
	    int itx,ity,itz;
	    FaceIterator fi=m.face.begin();
	    for (int i=0; i < faced ; i++) {
	      indicesf[fi] = i;
	      itx = it(0,i);
	      ity = it(1,i);
	      itz = it(2,i);
	      (*fi).V(0)=ivp[itx];
	      (*fi).V(1)=ivp[ity];
	      (*fi).V(2)=ivp[itz];
	      ++fi;
	    }
	    return 0;
	  } else {
	    return 1;
	  }
	} else {
	  return -1;
	}
      }
    };
}
