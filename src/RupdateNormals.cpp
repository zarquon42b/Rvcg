#include "typedef.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include <vcg/complex/algorithms/pointcloud_normal.h>

//using namespace vcg;
using namespace tri;
using namespace Rcpp;
//using namespace std;


RcppExport SEXP RupdateNormals(SEXP vb_, SEXP it_, SEXP type_, SEXP pointcloud_, SEXP silent_)
{
  try {
    // declare Mesh and helper variables
    int select = Rcpp::as<int>(type_);  
    Rcpp::IntegerVector pointcloud(pointcloud_);
    bool silent = as<bool>(silent_);
    MyMesh m;
    VertexIterator vi;
    FaceIterator fi;
    // allocate mesh and fill it
    int check = Rvcg::IOMesh<MyMesh>::RvcgReadR(m,vb_,it_);
    Rcpp::NumericMatrix normals(3,m.vn);
    if (check < 0) {
      ::Rf_error("mesh has no faces and/or no vertices");
      return wrap(1);
    } else if (check == 1) {
      if (!silent)
	Rprintf("%s\n","Info: mesh has no faces normals for point clouds are computed");
      PointCloudNormal<MyMesh>::Param p;
      p.fittingAdjNum = pointcloud[0];
      p.smoothingIterNum = pointcloud[1];
      p.viewPoint = Point3f(0,0,0);
      p.useViewPoint = false;
      PointCloudNormal<MyMesh>::Compute(m,p);
    }  else {
      // update normals
      if (select == 0) {
	tri::UpdateNormal<MyMesh>::PerVertex(m);
      } else {
	tri::UpdateNormal<MyMesh>::PerVertexAngleWeighted(m);
      }
      tri::UpdateNormal<MyMesh>::NormalizePerVertex(m);
     
      //write back
     
    }
    vi=m.vert.begin();
    SimpleTempData<MyMesh::VertContainer,int> indiceout(m.vert);
    for (int i=0;  i < m.vn; i++) {
      if( ! vi->IsD() )	{
	normals(0,i) = (*vi).N()[0];
	normals(1,i) = (*vi).N()[1];
	normals(2,i) = (*vi).N()[2];
      }
      ++vi;
    }
    
    return Rcpp::wrap(normals);
  
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
  
}
 

    
