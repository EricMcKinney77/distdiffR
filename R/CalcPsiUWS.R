#' The Psi UWS Statistic
#'
#' This statistic computes the uniformly weighted squared (UWS) differences between
#' the empirical cumulative distribution functions for the two samples.
#' For more information, see \insertCite{mckinney2022extensions;textual}{distdiffR}
#' and \insertCite{mckinney2021extensions;textual}{distdiffR}.
#'
#' @param data A two column matrix of the bivariate pooled samples
#' @param subjects A numerical vector of sample labels (use either 1 or 2)
#'
#' @return The Psi UWS statistic
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
#' CalcPsiUWS(as.matrix(pooled_data), sample_labels)
CalcPsiUWS <- function(data, subjects) {
  data1 <- data[subjects == 1, ]
  data2 <- data[subjects == 2, ]
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  sum((bcdf(data1, data1) - bcdf(data2, data1))^2) +
    sum((bcdf(data2, data2) - bcdf(data1, data2))^2)
}
