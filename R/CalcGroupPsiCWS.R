#' The Psi CWS statistic for aggregated group data
#'
#' This statistic computes the complementary weighted squared (CWS) differences
#' between the averaged subject empirical cumulative distribution functions for
#' the two samples. For more information, see
#' \insertCite{mckinney2022extensions;textual}{distdiffR} and
#' \insertCite{mckinney2021extensions;textual}{distdiffR}.
#'
#' @param data A two column matrix of the bivariate pooled samples
#' @param groups A numeric vector of sample (or group) labels (use either 1 or 2)
#' @param subjects A numeric vector of subject labels
#'
#' @return the Psi CWS statistic for aggregated group data
#' @export
#' @references
#' \insertRef{mckinney2022extensions}{distdiffR}
#'
#' \insertRef{mckinney2021extensions}{distdiffR}
#' @examples
#' # Randomly assign all three species to two samples
#' data(iris)
#' iris$Species <- rep(1:3, each = 50) # Species will serve as the subject label
#' irisPermuted <- iris[sample.int(nrow(iris)), ]
#' sample1 <- as.matrix(irisPermuted[1:75, c(1:2, 5)])
#' sample2 <- as.matrix(irisPermuted[76:150, c(1:2, 5)])
#' pooled_data <- rbind(cbind(sample1, 1), cbind(sample2, 2))
#'
#' CalcGroupPsiCWS(pooled_data[, 1:2], pooled_data[, 4], pooled_data[, 3])
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
