#ifndef __CHECK_LIST_NAMES__
#define __CHECK_LIST_NAMES__

#include <RcppArmadillo.h>
using Rcpp::List;

std::vector<bool> checkListNames(List mylist, Rcpp::CharacterVector mychar);
#endif
