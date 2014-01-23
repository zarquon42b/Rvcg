#include "typedef.h"
#include <wrap/ply/plylib.h>
#include <vcg/container/simple_temporary_data.h>
#include <wrap/io_trimesh/import.h>
#include <string.h>
#include "RvcgIO.h" 
#include <Rcpp.h>  
  
extern "C" {

  void RPlyRead(char **filename, double *vb ,int *dim, int *it, int *dimit, double *normals, int *getNorm, int *updNorm, double *quality,int *col, int *colvec, int *clean,int *fail)
  {
    int i;
    MyMesh m;
    // section read from input
    int faced = *dimit;
    //char file = **filename;
    char file[256];
    strcpy(file, *filename);
    int importNorm = *getNorm;
    int updateNorm = *updNorm;
    //load file
    int err2 = tri::io::ImporterPLY<MyMesh>::Open(m,file);
    if(err2) {
      Rprintf("Error in reading %s: '%s'\n",file,tri::io::Importer<MyMesh>::ErrorMsg(err2));
      //exit(-1);  
    }
    //Rprintf("%i",err2);
    if (err2 == 0)
      {
	if (updateNorm == 1) {
	  tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
	  tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
	}
	if (*clean == 1) {
	  int dup = tri::Clean<MyMesh>::RemoveDuplicateVertex(m);
	  int unref =  tri::Clean<MyMesh>::RemoveUnreferencedVertex(m);
	  Rprintf("Removed %i duplicate and %i unreferenced vertices\n",dup,unref);
	}
	vcg::tri::Allocator<MyMesh>::CompactVertexVector(m);
	vcg::tri::Allocator<MyMesh>::CompactFaceVector(m);
	
	if (m.vn > *dim || m.fn > *dimit) {
	  *fail = 1;
	  *dim=m.vn;
	  *dimit=m.fn;
	} else {
	  *dim=m.vn;
	  *dimit=m.fn;
	  
	  //--------------------------------------------------------------------------------------//
	  //                                   WRITE BACK
	  // Create meshes,
	  // Update the bounding box and initialize max search distance
	  // Remove duplicates and update mesh properties
	  //--------------------------------------------------------------------------------------//
	  SimpleTempData<MyMesh::VertContainer,int> indices(m.vert);
	  
	  //VertexPointer ivp[d];
	  if (m.vn > 0) {
	    VertexIterator vi=m.vert.begin();
	    for (i=0;  i < m.vn; i++) {
	      indices[vi] = i;//important: updates vertex indices
	      //	ivp[i]=&*vi;
	      vb[i*3] = (*vi).P()[0];
	      vb[i*3+1] = (*vi).P()[1];
	      vb[i*3+2] = (*vi).P()[2];
	      if (importNorm == 1) {
		normals[i*3] = (*vi).N()[0];
		normals[i*3+1] = (*vi).N()[1];
		normals[i*3+2] = (*vi).N()[2];
	      }
	      if (*col == 1) {
		colvec[i*3] = (*vi).C()[0];
		colvec[i*3+1] = (*vi).C()[1];
		colvec[i*3+2] = (*vi).C()[2];
	      }
	      ++vi;
	    }
	  }
	  
	  FacePointer fp;
	  int vv[3];
	  *dim = m.vn;
	  faced=m.fn;
	  
	  FaceIterator fi=m.face.begin();
	  
	  if (m.fn > 0) {
	    for (i=0; i < faced;i++) {
	      fp=&(*fi);
	      if( ! fp->IsD() ) {
		      vv[0]=indices[fp->cV(0)];
		      vv[1]=indices[fp->cV(1)];
		      vv[2]=indices[fp->cV(2)];
		      it[i*3]=vv[0];
		      it[i*3+1]=vv[1];
		      it[i*3+2]=vv[2];
		      ++fi;
	      }
	    }
	  }
	}
      }
  }
}

using namespace Rcpp;
RcppExport SEXP RPlyWrite(SEXP vb_, SEXP it_, SEXP binary_, SEXP addNormals_, SEXP filename_, SEXP colvec_, SEXP hasCol_)
{ 
  MyMesh m;
   //set up parameters 
  bool binary = Rcpp::as<bool>(binary_);
  bool addNormals = Rcpp::as<bool>(addNormals_);
  bool hasCol =  Rcpp::as<bool>(hasCol_);
  std::string str = Rcpp::as<std::string>(filename_);
  const char *filename = str.c_str();
  //char *filename[256] = strcpy(cstr
  //strcpy(filename1, filename);
  //allocate mesh and fill it
  Rcpp::IntegerMatrix colvec(colvec_);
  Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
  int mask0 = 0;
  
  if (addNormals) {
    tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
    mask0 = mask0 + tri::io::Mask::IOM_VERTNORMAL;
  }
  if (hasCol) {
    mask0 =mask0+ tri::io::Mask::IOM_VERTCOLOR;
    if (m.vn > 0) {
      int i;
      VertexIterator vi=m.vert.begin();
      for (i=0;  i < m.vn; i++) {
	(*vi).C()[0] = colvec(0, i);
	(*vi).C()[1] = colvec(1, i);
	(*vi).C()[2] = colvec(2, i);
	++vi;
      }
     
    }
    
    
  }
  
    tri::io::ExporterPLY<MyMesh>::Save(m, filename, mask0, binary);
  return Rcpp::wrap(0);
}

