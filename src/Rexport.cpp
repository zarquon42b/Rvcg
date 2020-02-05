#include "typedef.h"
#include <wrap/ply/plylib.h>
#include <vcg/container/simple_temporary_data.h>
#include <wrap/io_trimesh/import.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>
#include <wrap/io_trimesh/export_vrml.h>


#include <string.h>
#include "RvcgIO.h" 
#include <RcppArmadillo.h>  
  

using namespace Rcpp;
RcppExport SEXP RMeshWrite(SEXP mesh_, SEXP binary_, SEXP addNormals_, SEXP filename_, SEXP colvec_, SEXP hasCol_, SEXP writeNormals_, SEXP type_)
{ 
  try {
    MyMesh m;
    //set up parameters 
    List mesh(mesh_);
    bool binary = Rcpp::as<bool>(binary_);
    bool addNormals = Rcpp::as<bool>(addNormals_);
    bool hasCol =  Rcpp::as<bool>(hasCol_);
    bool writeNormals =  Rcpp::as<bool>(writeNormals_);
    int type = as<int>(type_);
    std::string str = Rcpp::as<std::string>(filename_);
    bool hasFaces = true;
    const char *filename = str.c_str();
    //char *filename[256] = strcpy(cstr
    //strcpy(filename1, filename);
    //allocate mesh and fill it
    Rcpp::IntegerMatrix colvec(colvec_);
    if (addNormals || writeNormals) {
      //m.vert.EnableNormal();
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
  
    Rvcg::IOMesh<MyMesh>::RvcgReadR(m,mesh["vb"],mesh["it"],mesh["normals"]);
    int mask0 = 0;
    
    if (addNormals) {
      if (hasFaces) {
	tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
      } else {
	Rcpp::IntegerVector pointcloud= IntegerVector::create(10,0);
	PointCloudNormal<MyMesh>::Param p;
	p.fittingAdjNum = pointcloud[0];
	p.smoothingIterNum = pointcloud[1];
	p.viewPoint = Point3f(0,0,0);
	p.useViewPoint = false;
	PointCloudNormal<MyMesh>::Compute(m,p);
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
    if (type == 0)
      tri::io::ExporterPLY<MyMesh>::Save(m, filename, mask0, binary);
    if (type == 1)
      tri::io::ExporterOFF<MyMesh>::Save(m, filename, mask0);
    if (type == 2)
      tri::io::ExporterOBJ<MyMesh>::Save(m, filename, mask0);
    if (type == 3)
      tri::io::ExporterSTL<MyMesh>::Save(m, filename, binary,mask0);
    if (type == 4)
      tri::io::ExporterWRL<MyMesh>::Save(m, filename, mask0,0);
    
    return Rcpp::wrap(0);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

// using namespace Rcpp;
// RcppExport SEXP RSTLWrite(SEXP vb_, SEXP it_, SEXP binary_, SEXP filename_)
// { 
// try {
//   MyMesh m;
//   //set up parameters 
//   bool binary = Rcpp::as<bool>(binary_);
//   std::string str = Rcpp::as<std::string>(filename_);
//   const char *filename = str.c_str();
//   //allocate mesh and fill it
//   Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    
//   tri::io::ExporterSTL<MyMesh>::Save(m, filename, binary);
//   return Rcpp::wrap(0);
// } catch (std::exception& e) {
//     ::Rf_error( e.what());
//     return wrap(1);
//   } catch (...) {
//     ::Rf_error("unknown exception");
//  }
// }


