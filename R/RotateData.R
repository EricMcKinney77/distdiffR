#' Create rotated versions of the data
#'
#' This function creates rotated versions of the pooled samples
#'
#' @param data a two column matrix of the bivariate pooled samples
#' @param numRotations a non-negative integer specifying the number of rotations to be applied to the data within 360 degrees.
#'
#' @return a list of matrices containing the coordinates for each version of the rotated data (including the original data)
#' @export
RotateData <- function(data, numRotations) {
  rotDataList <- vector("list", length = numRotations)
  rotDataList[[1]] <- data
  radRot <- 2 * pi * 0:(numRotations - 1) / numRotations
  for (i in 2:numRotations) {
    radRoti <- radRot[i] # Rotation angle in radians
    sinr <- sin(radRoti)
    cosr <- cos(radRoti)
    rotMat <- matrix(c(cosr, -sinr, sinr, cosr), 2, 2)
    rotDataList[[i]] <- data %*% rotMat
  }
  rotDataList
}
