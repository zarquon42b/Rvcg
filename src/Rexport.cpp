#include "typedefImport.h"
#include <wrap/ply/plylib.h>
#include <vcg/container/simple_temporary_data.h>
#include <wrap/io_trimesh/import.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>

#include <string.h>
#include "RvcgIO.h" 
#include <RcppArmadillo.h>  
  
/* extern "C" {

   void RPlyRead(char **filename, double *vb ,int *dim, int *it, int *dimit, double *normals, int *getNorm, int *updNorm, double *quality,int *col, int *colvec, int *clean,int *fail)
   {
   int i;
   MyMeshImport m;
   // section read from input
   int faced = *dimit;
   //char file = **filename;
   char file[256];
   strcpy(file, *filename);
   int importNorm = *getNorm;
   int updateNorm = *updNorm;
   //load file
   int err2 = tri::io::ImporterPLY<MyMeshImport>::Open(m,file);
   if(err2) {
   Rprintf("Error in reading %s: '%s'\n",file,tri::io::Importer<MyMeshImport>::ErrorMsg(err2));
   //exit(-1);  
   }
   //Rprintf("%i",err2);
   if (err2 == 0)
   {
   if (updateNorm == 1) {
   tri::UpdateNormal<MyMeshImport>::PerVertexAngleWeighted(m);
   tri::UpdateNormal<MyMeshImport>::NormalizePerVertex(m);
   }
   if (*clean == 1) {
   int dup = tri::Clean<MyMeshImport>::RemoveDuplicateVertex(m);
   int unref =  tri::Clean<MyMeshImport>::RemoveUnreferencedVertex(m);
   Rprintf("Removed %i duplicate and %i unreferenced vertices\n",dup,unref);
   }
   vcg::tri::Allocator<MyMeshImport>::CompactVertexVector(m);
   vcg::tri::Allocator<MyMeshImport>::CompactFaceVector(m);
	
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
   SimpleTempData<MyMeshImport::VertContainer,int> indices(m.vert);
	  
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
*/
using namespace Rcpp;
RcppExport SEXP RPlyWrite(SEXP mesh_, SEXP binary_, SEXP addNormals_, SEXP filename_, SEXP colvec_, SEXP hasCol_, SEXP writeNormals_)
{ 
  try {
    MyMeshImport m;
    //set up parameters 
    List mesh(mesh_);
    bool binary = Rcpp::as<bool>(binary_);
    bool addNormals = Rcpp::as<bool>(addNormals_);
    bool hasCol =  Rcpp::as<bool>(hasCol_);
    bool writeNormals =  Rcpp::as<bool>(writeNormals_);
    std::string str = Rcpp::as<std::string>(filename_);
    bool hasFaces = true;
    const char *filename = str.c_str();
    //char *filename[256] = strcpy(cstr
    //strcpy(filename1, filename);
    //allocate mesh and fill it
    Rcpp::IntegerMatrix colvec(colvec_);
    if (addNormals || writeNormals) {
      m.vert.EnableNormal();
    }
    Rcpp::CharacterVector normname("normals");
    Rcpp::CharacterVector nam = mesh.names();
    Rcpp::IntegerVector ind(Rf_match(nam,normname,0));
    Rcpp::LogicalVector   log(ind);
    
    if (log[0] == 0) {
      mesh["normals"] = wrap(0);
      writeNormals = false;
    }
    if (!Rf_isMatrix(mesh["it"]))
      hasFaces = false;
  
    Rvcg::IOMesh<MyMeshImport>::RvcgReadR(m,mesh["vb"],mesh["it"],mesh["normals"]);
    int mask0 = 0;
    
    if (addNormals) {
      if (hasFaces) {
	tri::UpdateNormal<MyMeshImport>::PerVertexAngleWeighted(m);
      } else {
	Rcpp::IntegerVector pointcloud= IntegerVector::create(10,0);
	PointCloudNormal<MyMeshImport>::Param p;
	p.fittingAdjNum = pointcloud[0];
	p.smoothingIterNum = pointcloud[1];
	p.viewPoint = Point3f(0,0,0);
	p.useViewPoint = false;
	PointCloudNormal<MyMeshImport>::Compute(m,p);
      }
	}
    if ( addNormals || writeNormals)
	 mask0 = mask0 + tri::io::Mask::IOM_VERTNORMAL;
    
    if (hasCol) {
      m.vert.EnableColor();
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
  
    tri::io::ExporterPLY<MyMeshImport>::Save(m, filename, mask0, binary);
    return Rcpp::wrap(0);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

using namespace Rcpp;
RcppExport SEXP RSTLWrite(SEXP vb_, SEXP it_, SEXP binary_, SEXP filename_)
{ 
try {
  MyMeshImport m;
  //set up parameters 
  bool binary = Rcpp::as<bool>(binary_);
  std::string str = Rcpp::as<std::string>(filename_);
  const char *filename = str.c_str();
  //allocate mesh and fill it
  Rvcg::IOMesh<MyMeshImport>::RvcgReadR(m,vb_,it_);
    
  tri::io::ExporterSTL<MyMeshImport>::Save(m, filename, binary);
  return Rcpp::wrap(0);
} catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
 }
}
