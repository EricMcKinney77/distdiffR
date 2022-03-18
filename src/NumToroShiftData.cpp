#include <Rcpp.h>
using namespace Rcpp;

//' Apply a toroidal shift to the pooled samples using a number of points
//'
//' The `NumToroShiftData()` function produces a list of toroidal shifted
//' versions of the two-column input matrix. The number of toroidal shifts
//' is an integer passed to `numShifts`. The origins of the toroidal shifts are
//' randomly selected from the combined samples. The pooled `data` is assumed to
//' list all of the first sample of size `n1` before the second sample (of size
//' `n2`).
//'
//' @param data A two column matrix of the pooled samples
//' @param n1 An integer sample size for the first sample
//' @param n2 An integer sample size for the second sample
//' @param numShifts A numeric number of points to be used as toroidal shift origins
//' @return A list of toroidal shifted pooled sample matrices
//' @examples
//' data(iris)
//' sample1 <- as.matrix(iris[iris$Species == "setosa", 1:2])
//' sample2 <- as.matrix(iris[iris$Species == "virginica", 1:2])
//' pooled_data <- rbind(sample1, sample2)
//' n1 <- nrow(sample1)
//' n2 <- nrow(sample2)
//'
//' # Create a list of five toroidal shifts of the pooled data
//' output <- NumToroShiftData(pooled_data, n1, n2, 25)
//' summary(output)
//' @export
// [[Rcpp::export]]
List NumToroShiftData(NumericMatrix data, int n1, int n2, int numShifts) {
  NumericMatrix::Column X = data( _, 0);
  NumericMatrix::Column Y = data( _, 1);
  float dataWidth = max(X) - min(X);
  float dataHeight = max(Y) - min(Y);

  int n1plusn2 = n1 + n2;
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
