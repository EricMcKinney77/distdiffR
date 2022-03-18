#' The Psi DWA Statistic
#'
#' This statistic computes the double weighted absolute (DWA) differences between
#' the empirical cumulative distribution functions for the two samples.
#' For more information, see \insertCite{mckinney2022extensions;textual}{distdiffR}
#' and \insertCite{mckinney2021extensions;textual}{distdiffR}.
#'
#' @param data a two column matrix of the bivariate pooled samples
#' @param subjects a numerical vector of sample labels (use either 1 or 2)
#'
#' @return the Psi DWA statistic
#' @export
#' @references
#' \insertRef{mckinney2022extensions}{distdiffR}
#'
#' \insertRef{mckinney2021extensions}{distdiffR}
#' @examples
#' data(iris)
#' pooled_data <- iris[iris$Species %in% c("setosa", "virginica"), 1:2]
#' sample_labels <- rep(1:2, c(sum(iris$Species == "setosa"),
#'                             sum(iris$Species == "virginica")))
#'
#' CalcPsiDWA(as.matrix(pooled_data), sample_labels)
CalcPsiDWA <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n1 * sum(abs(bcdf(data1, data1) - bcdf(data2, data1))) +
      n2 * sum(abs(bcdf(data2, data2) - bcdf(data1, data2)))) / (n1 + n2)
}
