#include <checkListNames.h>

std::vector<bool> checkListNames(List mylist, Rcpp::CharacterVector mychar) {
  try {
    
    Rcpp::CharacterVector nam = mylist.names();
    Rcpp::IntegerVector ind(Rf_match(nam,mychar,0));
    Rcpp::LogicalVector   log(ind);
    std::vector<bool> out = Rcpp::as<std::vector<bool> >(log);
    return out;
  } catch (std::exception& e) {
    ::Rf_error( e.what());
  } catch (...) {
    ::Rf_error("unknown exception");
  }
}
