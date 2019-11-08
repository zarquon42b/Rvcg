#include <vector>
#include <limits>
#include <RcppArmadillo.h>
#include <checkListNames.h>
#include <vcg/complex/complex.h>
#include <vcg/complex/append.h>
#include <vcg/container/simple_temporary_data.h>
#include <vcg/space/point3.h>

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
    static int RvcgReadR(MeshType &m, SEXP vb_, SEXP it_= Rcpp::wrap(0), SEXP normals_ = Rcpp::wrap(0), bool zerobegin=true, bool readnormals = true, bool readfaces=true) {
      try {
	//insert vertices
	if (Rf_isMatrix(vb_) && VertexType::HasCoord() ) {
	  Rcpp::NumericMatrix vb(vb_);
	  int d =  vb.ncol();
	  vcg::tri::Allocator<MeshType>::AddVertices(m,d);
	  std::vector<VertexPointer> ivp;
	  ivp.resize(d);
	  vcg::SimpleTempData<typename MeshType::VertContainer, unsigned int> indices(m.vert);
	  //read vertices
	  for (int i=0; i < d; i++) {
	    VertexIterator vi = m.vert.begin()+i;
	    ivp[i]=&*vi;
	    (*vi).P() = CoordType(vb(0,i),vb(1,i),vb(2,i));
	  }
	  //insert vertex normals
	  if (Rf_isMatrix(normals_) && vcg::tri::HasPerVertexNormal(m) && readnormals) {
	    Rcpp::NumericMatrix normals(normals_);
	    if (normals.ncol() != d) {
	      Rprintf("number of normals is not equal to number of vertices");
	    } else {
	      vcg::SimpleTempData<typename MeshType::VertContainer, unsigned int> indices(m.vert);
	      
	      // #pragma omp parallel for schedule(static)
	      for (int i=0; i < d; i++) {
		VertexIterator vi = m.vert.begin()+i;
		ivp[i]=&*vi;
		(*vi).N() = CoordType(normals(0,i),normals(1,i),normals(2,i));
	      }
	    }
	  }
	  //process faces but check attributes and input first
	  if (Rf_isMatrix(it_) && FaceType::HasVertexRef() && readfaces) {
	    Rcpp::IntegerMatrix it(it_);
	    unsigned int faced = it.ncol();
	    vcg::tri::Allocator<MeshType>::AddFaces(m,faced);
	    vcg::SimpleTempData<typename MeshType::FaceContainer, unsigned int> indicesf(m.face);
	    
	    for (unsigned int i=0; i < faced ; i++) {
	      int subtract = 0;
	      if (!zerobegin)
		subtract=1;
	      FaceIterator fi=m.face.begin()+i;
	      indicesf[fi] = i;
	      for (int j = 0; j < 3; j++) 
		(*fi).V(j)=ivp[it(j,i)-subtract];
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
	vcg::SimpleTempData<typename MeshType::VertContainer,unsigned int> indices(m.vert);
	Rcpp::NumericMatrix vb(4, m.vn), normals(4, m.vn);
	std::fill(vb.begin(),vb.end(),1);
	std::fill(normals.begin(),normals.end(),1);
	Rcpp::IntegerMatrix itout(3, m.fn);

	for (int i=0;  i < m.vn; i++) {
	  VertexIterator vi=m.vert.begin()+i;
	  indices[vi] = i;//important: updates vertex indices
	  for (int j = 0; j < 3; j++) {
	    vb(j,i) = (*vi).P()[j];
	    if (exnormals) 
	      normals(j,i) = (*vi).N()[j];
	  }
	}
	
	for (int i=0; i < m.fn;i++) {
	  FacePointer fp;
	  FaceIterator fi=m.face.begin()+i;
	  fp=&(*fi);
	  if (fp) {
	    if( ! fp->IsD() ) {
	      if (fp->V(0) && fp->V(1) && fp->V(2)) {
		for (int j = 0; j < 3; j++) {
		  itout(j,i) = indices[fp->cV(j)]+1;
		}
	      }
	    }
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
    
    static int mesh3d2Rvcg(MeshType &m, SEXP mesh_,bool zerobegin=false,bool readnormals=true,bool readfaces=true) {
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
      int out = RvcgReadR(m , mesh["vb"],mesh["it"],mesh["normals"],zerobegin,readnormals,readfaces);
      return out;
    };
    static arma::mat GetVertsArma(MeshType &m) {
      arma::mat vb(m.vn,3); 

      for (int i = 0; i < m.vn; i++) {
	VertexIterator vi=m.vert.begin()+i;
	for (int j = 0; j < 3; j++) 
	  vb(i,j) = (*vi).P()[j];
      }
      return vb;
    };
    static void VertsArmaToMesh(MeshType &m, const arma::mat &coords) {
      unsigned int d =  coords.n_rows;
      vcg::tri::Allocator<MeshType>::AddVertices(m,d);
      std::vector<VertexPointer> ivp;
      
      for (unsigned int i=0; i < d; i++) {
	VertexIterator vi = m.vert.begin()+i;
	(*vi).P() = CoordType(coords(i,0),coords(i,1),coords(i,2));
      }
    };
  };
}
