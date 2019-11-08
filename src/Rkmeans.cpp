#include "pointcloud.h"
#include "typedef.h"
#include "RvcgIO.h"
#include "RcppArmadillo.h"

using namespace Rcpp;
using namespace arma;

List fastSubsetMeans(mat &x, uvec &inds, int k, int threads) {
  try {

    mat center(k,x.n_cols);
    vec checkempty(k);checkempty.fill(0);
  
    center.fill(0);
#pragma omp parallel for schedule(static) num_threads(threads)
    for (int i =0; i < k;i++) {
      uvec tmpinds = find(inds ==i);
      mat tmpmat = x.rows(tmpinds);
      rowvec tmpresult(x.n_cols);tmpresult.fill(0);
      if (tmpinds.size() == 0)
	checkempty(i) = 1;
      for (unsigned int j = 0; j < tmpinds.size();j++)
	tmpresult += tmpmat.row(j);
      tmpresult /= tmpinds.size();
      center.row(i) = tmpresult;
    }
    List out = List::create(Named("centers")=center,
			    Named("checkempty")=checkempty);
				
    return out;
  }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}


RcppExport SEXP Rkmeans(SEXP mesh_, SEXP k_, SEXP itermax_, SEXP threads_) {
  try {
  PcMesh mesh;
  int k = as<int>(k_);
  int itermax = as<int>(itermax_);
  int threads = as<int>(threads_);
  int nofPointsPerCell = 16;
  int maxDepth = 64;
  Rvcg::IOMesh<PcMesh>::mesh3d2Rvcg(mesh,mesh_,false,false);
  arma::mat coords = Rvcg::IOMesh<PcMesh>::GetVertsArma(mesh);
  int npts = coords.n_rows;
  uvec samplevec(npts);
  uvec subset(k);
  //initialize index vectors
  for (int i =0; i < npts;i++) {
    samplevec[i] = i;
    if (i < k)
      subset[i] = i;
  }
  
  uvec shufflesample = shuffle(samplevec);
  uvec subset2 = shufflesample(subset);
  mat centers = coords.rows(subset2);
  int count = 0;
  double centercheck = 1e12;
  //set up kdtree search
  VertexConstDataWrapper<PcMesh> ww(mesh);
  uvec clostinds = samplevec;
  
  while (count < itermax && centercheck != 0) {
    uvec clost_old = clostinds;
    PcMesh centermesh;
    Rvcg::IOMesh<PcMesh>::VertsArmaToMesh(centermesh,centers);
    VertexConstDataWrapper<PcMesh> ww(centermesh);
    KdTree<float> kdtree(ww, nofPointsPerCell, maxDepth);
    KdTree<float>::PriorityQueue queue;

#pragma omp parallel for schedule(static) firstprivate(queue, kdtree) num_threads(threads)
    for (int i = 0; i < mesh.vn; i++) {
      kdtree.doQueryK(mesh.vert[i].cP(), 1, queue);
      int neighbours = queue.getNofElements();
      clostinds[i] = queue.getIndex(0);
    }
    
    List subsetmeans = fastSubsetMeans(coords,clostinds,k,threads);
    centers = as<mat>(subsetmeans["centers"]);
    vec checkempty = as<vec>(subsetmeans["checkempty"]);
    //check if there are empty clusters, reshuffle and start again
    if (arma::sum(checkempty > 0)) {
	shufflesample = shuffle(samplevec);
	clostinds = shufflesample;
	subset2 = shufflesample(subset);
	centers = coords.rows(subset2);
	count = 0;
      }
    centercheck = arma::max(arma::abs(clostinds-clost_old));
    count += 1;
    
  }
  //clostinds -= 1;
  return List::create(Named("centers") = centers,
		      Named("class") = clostinds);
   }  catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
