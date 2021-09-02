#include <Rcpp.h>
using namespace Rcpp;

//' Apply a toroidal shift to the pooled samples using a proportion of points
//'
//' @param data A two column matrix of the pooled samples
//' @param n1 An integer sample size for the first sample
//' @param n2 An integer sample size for the second sample
//' @param propPnts A numeric proportion of points to be used as toroidal shift origins
//' @return A two column matrix of the shifted pooled samples
//' @export
// [[Rcpp::export]]
List PropToroShiftData(NumericMatrix data, int n1, int n2, float propPnts = 1) {
  NumericMatrix::Column X = data( _, 0);
  NumericMatrix::Column Y = data( _, 1);
  float dataWidth = max(X) - min(X);
  float dataHeight = max(Y) - min(Y);

  int n1plusn2 = n1 + n2;
  int numShifts = round(propPnts * n1plusn2);
  List shiftDataList = List::create();
  IntegerVector n1n2vec = seq(0, n1plusn2 - 1);
  IntegerVector sampledRows = sample(n1n2vec, numShifts);
  for (int x = 0; x < numShifts; ++x) {
    int sampledRow = sampledRows[x];
    NumericVector toroCntrPnt = data(sampledRow, _);
    NumericVector shiftedX = ifelse(X < toroCntrPnt[0], X + dataWidth, X);
    NumericVector shiftedY = ifelse(Y < toroCntrPnt[1], Y + dataHeight, Y);
    shiftDataList.push_back(cbind(shiftedX, shiftedY));
  }
  return shiftDataList;
}
