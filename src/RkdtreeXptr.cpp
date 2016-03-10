#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <RcppArmadillo.h>
#include <RvcgKD.h>

RcppExport SEXP createKDtree(SEXP target_, SEXP nofPointsPerCell_, SEXP maxDepth_) {
  Rcpp::XPtr< MyMesh > target(new MyMesh,true);
  Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(*target,target_,false,false);
  unsigned int nofPointsPerCell = as<unsigned int >(nofPointsPerCell_);
  unsigned int maxDepth = as<unsigned int >(maxDepth_);
  VertexConstDataWrapper<MyMesh> ww(*target);
  KdTree<float> tree(ww, nofPointsPerCell, maxDepth);

  Rcpp::XPtr< KdTree<float> > Rtree(new KdTree<float>(ww, nofPointsPerCell, maxDepth),true);
  
  
  return List::create(Named("kdtree") =Rtree,
		      Named("target") = target);

}

List searchKDtree(KdTree<float> kdtree, MyMesh &target, MyMesh &query, int k, int threads){
  try {
    KdTree<float>::PriorityQueue queue;
    typedef pair<float,int> mypair;
    IntegerMatrix result(query.vn,k);
    NumericMatrix distance(query.vn,k);
    std::fill(result.begin(), result.end(),-1);
#pragma omp parallel for firstprivate(queue, kdtree) schedule(static) num_threads(threads)
    for (int i = 0; i < query.vn; i++) {
      //tree.doQueryK(query.vert[i].cP());
      kdtree.doQueryK(query.vert[i].cP(), k, queue);
      //int neighbours = tree.getNofFoundNeighbors();
      int neighbours = queue.getNofElements();
      vector<mypair> sortit;
      for (int j=0; j < neighbours; j++) {      
	int neightId = queue.getIndex(j);
	float dist = Distance(query.vert[i].cP(),target.vert[neightId].cP());
	sortit.push_back(mypair(dist, neightId));
      }

      sort(sortit.begin(),sortit.end());
      for (int j = 0; j < neighbours; j++){
	result(i,j) = sortit[j].second;
	distance(i,j) = sortit[j].first;
      }
    }
    return List::create(Named("index") = result,
			Named("distance") = distance)
      ;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RsearchKDtree(SEXP kdtree_,SEXP target_, SEXP query_, SEXP k_, SEXP threads_) {
  try {
    XPtr< KdTree<float> > kdtree(kdtree_);
    XPtr< MyMesh > target(target_);
    MyMesh query;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(query,query_,false,false,false);
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    List out = searchKDtree(*kdtree, *target, query, k, threads);
    return out;
    // start searching
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}

RcppExport SEXP RsearchKDtreeForClosestPoints(SEXP kdtree_,SEXP bary_, SEXP targetmesh_, SEXP query_, SEXP k_,SEXP sign_,SEXP borderchk_,SEXP barycentric_, SEXP angdev_=wrap(0), SEXP wnorm_=wrap(true), SEXP facenormals_=wrap(false), SEXP threads_=wrap(1)) {
  try {
    XPtr< KdTree<float> > kdtree(kdtree_);
    XPtr< MyMesh > bary(bary_);
    XPtr< MyMesh > target(targetmesh_);
    MyMesh query;
    Rvcg::IOMesh<MyMesh>::mesh3d2Rvcg(query,query_,false,true,false);
    int k = as<int>(k_);
    int threads = as<int>(threads_);
    bool barycentric = as<bool>(barycentric_);
    bool borderchk = as<bool>(borderchk_);
    bool wnorm = as<bool>(wnorm_);
    double angdev = as<double>(angdev_);
    bool sign = as<bool>(sign_);
    bool facenormals = as<bool>(facenormals_);
    (*target).face.EnableNormal();
    if (angdev > 0) {
      tri::UpdateNormal<MyMesh>::PerVertexNormalized(query);
    }
    tri::UpdateNormal<MyMesh>::PerFaceNormalized(*target);
    tri::UpdateNormal<MyMesh>::PerVertexNormalized(*target);
    if (borderchk) { //update Border flags
      tri::UpdateFlags<MyMesh>::FaceBorderFromNone(*target);
      tri::UpdateSelection<MyMesh>::FaceFromBorderFlag(*target);
    }
    List indices = searchKDtree(*kdtree,*bary,query,k,threads);
    arma::imat closest_indices = as<arma::imat>(indices["index"]);

    List out = Rvcg::KDtree<MyMesh,MyMesh>::clostKD(*target, query, closest_indices, k,angdev, facenormals,sign, wnorm,borderchk,barycentric ,threads);
    return out;
    
  } catch (std::exception& e) {
    ::Rf_error( e.what());
    return wrap(1);
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
  


