// Author: Stefan Schlager
// Date: 15 September 2010

#include "typedef.h"
#include "RvcgIO.h" 
#include <Rcpp.h>

using namespace Rcpp;

RcppExport SEXP RvertAdj(SEXP _vn, SEXP _it)
{
  Rcpp::IntegerMatrix it(_it);
  int vn = Rcpp::as<int>(_vn);
  int faced = it.ncol();
  int i, j;
  Rcpp::List outlist(vn);

  for (i=0; i < vn; i++) {
    std::vector<int> tmp;
    for (j=0; j < faced; j++) {
      if (it(0,j) == i || it(1,j) ==i || it(2,j) == i)
	tmp.push_back(j+1);
    }
    outlist[i] = tmp;
  }
  
  
  return outlist;
}
  
   


