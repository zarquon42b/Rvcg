#include "typedef.h"
#include "RvcgIO.h"
#include <Rcpp.h>
#include<vcg/complex/algorithms/create/platonic.h>

using namespace tri;
using namespace Rcpp;

RcppExport SEXP Rplatonic(SEXP subdiv_ = wrap(3)) {
  int subdiv = as<int>(subdiv_);
  MyMesh m;
  Sphere(m,subdiv);
  List out = Rvcg::IOMesh<MyMesh>::RvcgToR(m,false);
  return out;
}
