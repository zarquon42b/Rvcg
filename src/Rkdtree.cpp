#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include <RvcgKD.h>
#include <Rconfig.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rkdtree(SEXP vb0_, SEXP vb1_, SEXP k_ ,SEXP nofP_= wrap(16),SEXP mDepth_= wrap(64),SEXP threads_=wrap(1)) {
  try {
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    unsigned int nofP = as<unsigned int >(nofP_);
    unsigned int mDepth = as<unsigned int >(mDepth_);
    typedef pair<float,int> mypair;
    PcMesh target, query;
    Rvcg::IOMesh<PcMesh>::RvcgReadR(target, vb0_);  
    Rvcg::IOMesh<PcMesh>::RvcgReadR(query, vb1_);
 
    List out = Rvcg::KDtree< PcMesh, PcMesh >::KDtreeIO(target, query, k,nofP, mDepth,threads);
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
  
}

RcppExport SEXP RclosestKD(SEXP target_, SEXP query_, SEXP k_, SEXP sign_, SEXP smooth_, SEXP barycentric_, SEXP borderchk_, SEXP nofP_= wrap(16),SEXP mDepth_= wrap(64),SEXP angdev_=wrap(0), SEXP wnorm_=wrap(true), SEXP facenormals_ = wrap(true), SEXP threads_=wrap(1)) {
  try {
    bool smooth = as<bool>(smooth_);
    bool barycentric = as<bool>(barycentric_);
    bool borderchk = as<bool>(borderchk_);
    bool wnorm = as<bool>(wnorm_);
    bool facenormals = as<bool>(facenormals_);
    unsigned int nofP = as<unsigned int >(nofP_);
    unsigned int mDepth = as<unsigned int >(mDepth_);
    int threads = as<int>(threads_);
    int k = as<int>(k_);
    bool sign = as<bool>(sign_);
    MyMesh target;
    PcMesh  bary;
    MyMesh query;
    int checkit = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(target,target_);
    double angdev = as<double>(angdev_);
    target.face.EnableNormal();
    checkit = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(query, query_);
    if (angdev > 0) {
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(query);
    }
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(target);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(target);
    if (smooth) {
      tri::Smooth<MyMesh>::VertexNormalLaplacian(target,2,false);
      tri::UpdateNormal<MyMesh>::NormalizePerVertex(target);
    }
    if (borderchk) { //update Border flags
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(target);
      tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(target);
    }
    Rvcg::KDtree< MyMesh, PcMesh >::getBary(target, bary);
    List indices = Rvcg::KDtree< PcMesh, MyMesh >::KDtreeIO(bary, query, k,nofP, mDepth,threads);
    arma::imat closest_indices = as<arma::imat>(indices["index"]);
    List out = Rvcg::KDtree<MyMesh,MyMesh>::clostKD(target, query, closest_indices, k,angdev, facenormals,sign, wnorm,borderchk,barycentric ,threads);
    return out;
      
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP Rbarycenter(SEXP mesh_) {
  try {
    MyMesh m;
    int checkit = Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(m,mesh_);
    MyMesh out;
    Rvcg::KDtree< MyMesh, MyMesh >::getBary(m, out);
    Rcpp::NumericMatrix barycoord(3,out.vn);
    for (int i = 0; i < out.vn; i++) {
      MyMesh::CoordType tmp;
      tmp = out.vert[i].cP();
      barycoord(0,i) = tmp[0];
      barycoord(1,i) = tmp[1];
      barycoord(2,i) = tmp[2];
    }
    return wrap(barycoord);
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
