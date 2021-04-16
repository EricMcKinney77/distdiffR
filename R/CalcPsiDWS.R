#' The Psi_DWS Statistic
#'
#' This statistic computes the double weighted squared (DWS) differences between
#' the empirical cumulative distribution functions for the two samples.
#'
#' @param data a two column matrix of the bivariate pooled samples
#' @param subjects a numerical vector of sample labels (use either 1 or 2)
#'
#' @return the Psi_DWS statistic
#' @export
#'
#' @examples
CalcPsiDWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n1 * sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
      n2 * sum((bcdf(data2, data2) - bcdf(data1, data2))^2)) / (n1 + n2)
}
