#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP threshold(SEXP array_,  SEXP lowerbound_,  SEXP upperbound_){
//, SEXP lowerbound_, SEXP upperbound_) {
  IntegerVector array(array_);
  int lb = as<int>(lowerbound_);
  int ub = as<int>(upperbound_);
  int arlen = array.size();
  for (int i =0; i < arlen;i++) {
    if (array[i] > lb && array[i] < ub) {
      array[i] = 1;
    } else {
      array[i] = 0;
    }
  }
  return wrap(array);
}

