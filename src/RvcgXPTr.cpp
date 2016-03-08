#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>


RcppExport SEXP RmeshXPtr(SEXP mesh_) {
   Rcpp::XPtr< MyMesh > mesh(new MyMesh,true);
   Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(*mesh,mesh_);
   return mesh;
}
