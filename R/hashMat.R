#' Assigns a hash value to a two-column matrix.
#'
#' The `hashMat()` function assigns hash numbers to the two-column sample
#' matrices for the purpose of providing identical test results regardless of
#' the order in which the input data is passed to the `distdiffr()` or
#' `grouped_distdiffr()` functions.
#'
#' @param mat A two column matrix of bivariate observations.
#'
#' @return A numeric hash value
#' @export
hashMat <- function(mat) {
  # Normalize the data
  matMins <- apply(mat, 2, min)
  matRanges <- as.vector(diff(apply(mat, 2, range)))
  mat <- scale(mat, center = matMins, scale = matRanges)

  # Store the dimensions of the matrix
  dims <- dim(mat)

  # Construct and return the hash value
  hash <- 2
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      hash <- hash / 3 + mat[i, j]
    }
  }
  hash
}
