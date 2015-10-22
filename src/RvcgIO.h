#include <vector>
#include <limits>
#include <Rcpp.h>
#include <checkListNames.h>
#include <vcg/complex/complex.h>
//#include <vcg/complex/allocate.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/space/point3.h>
// #include <Rconfig.h>
// #ifdef SUPPORT_OPENMP
// #include <omp.h>
// #endif
using Rcpp::List;

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
    static int RvcgReadR(MeshType &m, SEXP vb_, SEXP it_= Rcpp::wrap(0), SEXP normals_ = Rcpp::wrap(0)) {
      try {
	//insert vertices
	if (Rf_isMatrix(vb_) && VertexType::HasCoord() ) {
	  Rcpp::NumericMatrix vb(vb_);
	  size_t d =  vb.ncol();
	  ScalarType x,y,z;	 
	  vcg::tri::Allocator<MeshType>::AddVertices(m,d);
	  std::vector<VertexPointer> ivp;
	  ivp.resize(d);
	  vcg::SimpleTempData<typename MeshType::VertContainer, int> indices(m.vert);
	  //VertexIterator vi = m.vert.begin();
	  // #pragma omp parallel for schedule(static)
	  for (size_t i=0; i < d; i++) {
	    VertexIterator vi = m.vert.begin()+i;
	    ivp[i]=&*vi;
	    x = vb(0,i);
	    y = vb(1,i);
	    z = vb(2,i);
	    (*vi).P() = CoordType(x,y,z);
	  }
	  //insert vertex normals
	  if (Rf_isMatrix(normals_) && vcg::tri::HasPerVertexNormal(m)) {
	    Rcpp::NumericMatrix normals(normals_);
	    if (normals.ncol() != d) {
	      Rprintf("number of normals is not equal to number of vertices");
	    } else {
	      vcg::SimpleTempData<typename MeshType::VertContainer, int> indices(m.vert);
	      
	      // #pragma omp parallel for schedule(static)
	      for (size_t i=0; i < d; i++) {
		VertexIterator vi = m.vert.begin()+i;
		ivp[i]=&*vi;
		x = normals(0,i);
		y = normals(1,i);
		z = normals(2,i);
		(*vi).N() = CoordType(x,y,z);
		
	      }
	    }
	  }
	  //process faces but check attributes and input first
	  if (Rf_isMatrix(it_) && FaceType::HasVertexRef()) {
	    Rcpp::IntegerMatrix it(it_);
	    int faced = it.ncol();
	    vcg::tri::Allocator<MeshType>::AddFaces(m,faced);
	    vcg::SimpleTempData<typename MeshType::FaceContainer, int> indicesf(m.face);
	    
	    // #pragma omp parallel for schedule(static)
	    for (size_t i=0; i < faced ; i++) {
	      int itx,ity,itz;
	      FaceIterator fi=m.face.begin()+i;
	      indicesf[fi] = i;
	      itx = it(0,i);
	      ity = it(1,i);
	      itz = it(2,i);
	      (*fi).V(0)=ivp[itx];
	      (*fi).V(1)=ivp[ity];
	      (*fi).V(2)=ivp[itz];
	    }
	    return 0;
	  } else {
	    return 1;
	  }
	} else {
	  return -1;
	}
      } catch (std::exception& e) {
	::Rf_error( e.what());
      } catch (...) {
	::Rf_error("unknown exception");
      }
    };
   
    static Rcpp::List RvcgToR(MeshType &m, bool exnormals=false) {
      try {
	List out;
	vcg::SimpleTempData<typename MeshType::VertContainer,int> indices(m.vert);
	Rcpp::NumericMatrix vb(4, m.vn), normals(4, m.vn);
	std::fill(vb.begin(),vb.end(),1);
	std::fill(normals.begin(),normals.end(),1);
	Rcpp::IntegerMatrix itout(3, m.fn);
	// #pragma omp parallel for schedule(static)
	for (unsigned int i=0;  i < m.vn; i++) {
	  VertexIterator vi=m.vert.begin()+i;
	  indices[vi] = i;//important: updates vertex indices
	  vb(0,i) = (*vi).P()[0];
	  vb(1,i) = (*vi).P()[1];
	  vb(2,i) = (*vi).P()[2];
	  if (exnormals) {
	    normals(0,i) = (*vi).N()[0];
	    normals(1,i) = (*vi).N()[1];
	    normals(2,i) = (*vi).N()[2];
	  }
	}
	// #pragma omp parallel for schedule(static)
	for (size_t i=0; i < m.fn;i++) {
	  FacePointer fp;
	  FaceIterator fi=m.face.begin()+i;
	  fp=&(*fi);
	  if( ! fp->IsD() ) {
	    itout(0,i) = indices[fp->cV(0)]+1;
	    itout(1,i) = indices[fp->cV(1)]+1;
	    itout(2,i) = indices[fp->cV(2)]+1;
	  }
	}
	out["vb"] = vb;
	out["it"] = itout;
	if (exnormals)
	  out["normals"] = normals;
	out.attr("class") = "mesh3d";
	return out;
      
      } catch (std::exception& e) {
	::Rf_error( e.what());
      } catch (...) {
	::Rf_error("unknown exception");
      }
    };
    
    static int mesh3d2Rvcg(MeshType &m, SEXP mesh_) {
      List mesh(mesh_);
      Rcpp::CharacterVector mychar = Rcpp::CharacterVector::create("vb","it","normals");
      std::vector<bool> test = checkListNames(mesh,mychar);
      for (int i = 0; i < 3; i++) {
	if (!test[i]) {
	  std::string tmp = Rcpp::as<std::string>(mychar[i]);
	  mesh[tmp] = Rcpp::wrap(0);
	}
      }
      if (!test[0])
	::Rf_error("mesh has no vertices");
      int out = RvcgReadR(m , mesh["vb"],mesh["it"],mesh["normals"]);
      return out;
    };
  };
}
