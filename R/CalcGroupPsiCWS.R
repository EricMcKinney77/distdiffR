#' The Psi_CWS statistic for aggregated group data
#'
#' This statistic computes the complementary weighted squared (CWS) differences
#' between the averaged subject empirical cumulative distribution functions for
#' the two samples.
#'
#' @param data A two column matrix of the bivariate pooled samples
#' @param groups A numeric vector of sample (or group) labels (use either 1 or 2)
#' @param subjects A numeric vector of subject labels
#'
#' @return the Psi_CWS statistic
#' @export
CalcGroupPsiCWS <- function(data, groups, subjects) {
  data1 <- data[groups == 1, ]
  data2 <- data[groups == 2, ]
  subjects1 <- subjects[groups == 1]
  subjects2 <- subjects[groups == 2]
  subjIndices <- 1:max(subjects)
  data1bySubsLst <- lapply(subjIndices, function(x) data1[subjects1 == x, ])
  data2bySubsLst <- lapply(subjIndices, function(x) data2[subjects2 == x, ])
  bcdf11Lst <- lapply(subjIndices, function(x) bcdf(data1bySubsLst[[x]], data1))
  bcdf21Lst <- lapply(subjIndices, function(x) bcdf(data2bySubsLst[[x]], data1))
  bcdf22Lst <- lapply(subjIndices, function(x) bcdf(data2bySubsLst[[x]], data2))
  bcdf12Lst <- lapply(subjIndices, function(x) bcdf(data1bySubsLst[[x]], data2))
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  (n2 * sum((apply(do.call(cbind, bcdf11Lst), 1, mean) - apply(do.call(cbind, bcdf21Lst), 1, mean))^2) +
      n1 * sum((apply(do.call(cbind, bcdf22Lst), 1, mean) - apply(do.call(cbind, bcdf12Lst), 1, mean))^2)) / (n1 + n2)
}
