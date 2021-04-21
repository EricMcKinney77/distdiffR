#include <Rcpp.h>
using namespace Rcpp;

//' Construct and compute bivariate cumulative distribution function values
//'
//' @param data A two column matrix for constructing the empirical bivariate cumulative distribution function (EBCDF)
//' @param eval A two column matrix for input into the EBCDF
//' @return A numeric vector of output values from the EBCDF
//' @export
// [[Rcpp::export]]
NumericVector bcdf(NumericMatrix data, NumericMatrix eval) {
  NumericVector x(eval.nrow());
  int i, j;
  for(i = 0; i < eval.nrow(); ++i) {
    x(i) = 0;
    for(j = 0; j < data.nrow(); ++j) {
      if( data(j, 0) <= eval(i, 0) && data(j, 1) <= eval(i, 1) ) x(i) = x(i) + 1.0;
    }
    x(i) = x(i) / data.nrow();
  }
  return x;
}
