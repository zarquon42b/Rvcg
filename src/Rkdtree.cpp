#include "typedef.h"
#include "pointcloud.h"
#include "RvcgIO.h"
#include <Rcpp.h>

using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rkdtree(SEXP vb0_, SEXP vb1_, SEXP k_) {
  int k = as<int>(k_);
  typedef pair<float,int> mypair;
  PcMesh target, query;
  Rvcg::IOMesh<PcMesh>::RvcgReadR(target, vb0_);  
  Rvcg::IOMesh<PcMesh>::RvcgReadR(query, vb1_);
  VertexConstDataWrapper<PcMesh> ww(target);
  IntegerMatrix result(query.vn,k);
  std::fill(result.begin(), result.end(),-1);
  NumericMatrix distance(query.vn,k);
  KdTree<float> tree(ww);
  tree.setMaxNofNeighbors(k);
  for (int i = 0; i < query.vn; i++) {
    tree.doQueryK(query.vert[i].cP());
    int neighbours = tree.getNofFoundNeighbors();
    vector<mypair> sortit;
    for (int j=0; j < neighbours; j++) {      
      int neightId = tree.getNeighborId(j);
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
  }
