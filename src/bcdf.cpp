#include <Rcpp.h>
using namespace Rcpp;

//' Construct and evaluate a bivariate empirical cumulative distribution function
//'
//' Construct a bivariate empirical cumulative distribution function (BECDF)
//' using `data` and pass each of the `eval` points through the BECDF.
//'
//' @param data A two column matrix for constructing the BECDF
//' @param eval A two column matrix for input into the BECDF
//' @return A numeric vector of output values from the BECDF
//' @examples
//' data(iris)
//' sample1 <- as.matrix(iris[iris$Species == "virginica", 1:2])
//' sample2 <- as.matrix(iris[iris$Species == "versicolor", 1:2])
//'
//' bcdf(sample1, sample2)
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
