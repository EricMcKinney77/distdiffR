#' The Psi_UWA Statistic
#'
#' This statistic computes the uniformly weighted absolute (UWA) differences between
#' the empirical cumulative distribution functions for the two samples.
#'
#' @param data a two column matrix of the bivariate pooled samples
#' @param subjects a numerical vector of sample labels (use either 1 or 2)
#'
#' @return the Psi_UWA statistic
#' @export
CalcPsiUWA <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  sum(abs(bcdf(data1, data1) - bcdf(data2, data1))) +
    sum(abs(bcdf(data2, data2) - bcdf(data1, data2)))
}
